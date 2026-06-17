"""
Dialogs for displaying ESL-PSC analysis results.
"""
from __future__ import annotations
import os

import math
import numpy as np
# Use shared FASTA reader
from gui.core.fasta_io import read_fasta
import pandas as pd
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QBrush, QFontMetrics, QKeySequence, QGuiApplication, QShortcut
from PySide6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QSizePolicy, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QTreeWidget, QTreeWidgetItem, QMessageBox,
    QLabel, QHeaderView, QPushButton, QHBoxLayout,
    QComboBox, QDialogButtonBox, QFileDialog, QListWidget, QProgressDialog
)
from PySide6.QtSvgWidgets import QSvgWidget

from gui.ui.widgets.protein_map import ProteinMapWidget
from gui.ui.widgets.site_viewer import SiteViewer
from gui.ui.widgets.dialogs import PhenoThresholdDialog

# Keep references to open dialogs to prevent garbage collection
_open_dialogs = []


class NumericItem(QTableWidgetItem):
    """A QTableWidgetItem that sorts by numeric value while displaying formatted text."""
    def __init__(self, value: float, display_text: str):
        super().__init__(display_text)
        try:
            self._value = float(value)
        except Exception:
            self._value = float('nan')
        # Align numbers to the right for readability
        self.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

    def __lt__(self, other):
        if isinstance(other, NumericItem):
            # Handle NaNs by pushing them to the end
            if self._value != self._value:
                return False
            if other._value != other._value:
                return True
            return self._value < other._value
        # Fallback to default behavior
        try:
            return float(self.text()) < float(other.text())
        except Exception:
            return super().__lt__(other)

def _launch_site_viewer(
    gene: str,
    config,
    sites_path: str | None,
    parent=None,
    outgroup_species: list[str] | None = None,
    prefer_ccs_filter: bool = False,
) -> None:
    """Open the SiteViewer window for the given gene."""
    align_dir = getattr(config, "alignments_dir", "")
    if not align_dir:
        raise RuntimeError("Alignment directory not specified")
    align_path = os.path.join(align_dir, f"{gene}.fas")
    if not os.path.exists(align_path):
        raise FileNotFoundError(f"Alignment file not found: {align_path}")

    # Load alignment records using shared FASTA reader for reliability/consistency
    records = read_fasta(align_path)

    # Determine species groups with the following priority:
    # 1) If a response_dir with matrices is present, use the FIRST matrix file to define groups.
    #    This supports flipped-null and multimatrix runs where group membership per combo
    #    comes from the response matrices rather than the species groups file.
    # 2) Otherwise, fall back to the species groups file (first two non-empty lines only).
    conv: list[str] = []
    ctrl: list[str] = []

    # If parent dialog has an explicitly selected groups combo, respect it first.
    # This enables Site Counter to open SiteViewer with a user-chosen species grouping
    # derived from the species groups file.
    try:
        if parent is not None and hasattr(parent, '_selected_groups_combo'):
            sel = getattr(parent, '_selected_groups_combo')
            if isinstance(sel, (tuple, list)) and len(sel) == 2:
                conv = list(sel[0])
                ctrl = list(sel[1])
                used_response_dir = True  # treat as explicit selection
            else:
                conv = []
                ctrl = []
    except Exception:
        conv = []
        ctrl = []

    if not conv and not ctrl:
        try:
            pref = getattr(config, 'preferred_groups_combo', None)
            if pref and isinstance(pref, (tuple, list)) and len(pref) == 2:
                conv = list(pref[0])
                ctrl = list(pref[1])
        except Exception:
            pass

    # Prefer response_dir if available (and no explicit groups combo selected)
    response_dir = getattr(config, "response_dir", "") or ""
    if not response_dir and not conv and not ctrl:
        # Derive default response_dir name from species_groups_file if not explicitly set
        base = os.path.basename(getattr(config, "species_groups_file", "")).replace(".txt", "")
        if base and getattr(config, "output_dir", ""):
            response_dir = os.path.join(config.output_dir, f"{base}_response_matrices")

    used_response_dir = False
    resp_values: dict[str, float] = {}
    try:
        if (not conv and not ctrl) and response_dir and os.path.isdir(response_dir):
            files = sorted([f for f in os.listdir(response_dir) if f.endswith('.txt')])
            if files:
                # Determine which matrix to use: prefer a selection stored on the parent dialog,
                # otherwise default to the first file. Never prompt here.
                selected_name = files[0]
                try:
                    if parent is not None and hasattr(parent, '_selected_response_matrix'):
                        sel = getattr(parent, '_selected_response_matrix')
                        if sel in files:
                            selected_name = sel
                    elif getattr(config, 'preferred_response_matrix', '') in files:
                        selected_name = getattr(config, 'preferred_response_matrix')
                except Exception:
                    pass

                first_matrix = os.path.join(response_dir, selected_name)
                # Parse by line order: even index -> Convergent, odd index -> Control.
                # This holds for both binary and continuous response matrices.
                with open(first_matrix, encoding="utf-8", errors="ignore") as f:
                    for idx, raw in enumerate(f):
                        line = raw.strip()
                        if not line:
                            continue
                        parts = line.split()  # tab or whitespace separated
                        if len(parts) < 2:
                            continue
                        sp = parts[0]
                        # capture numeric phenotype value if present (binary or continuous)
                        try:
                            resp_values[sp] = float(parts[1])
                        except ValueError:
                            pass
                        if idx % 2 == 0:
                            conv.append(sp)
                        else:
                            ctrl.append(sp)
                used_response_dir = True
                if not getattr(config, 'preferred_response_matrix', ''):
                    try:
                        config.preferred_response_matrix = selected_name
                    except Exception:
                        pass
    except Exception:
        # Fall through to species groups parsing on error
        conv, ctrl, used_response_dir = [], [], False

    if not used_response_dir and not conv and not ctrl:
        # Fall back to species groups file: construct the FIRST combo by
        # taking the first species from each even-indexed (convergent) line
        # and each odd-indexed (control) line.
        groups_path = getattr(config, "species_groups_file", "")
        if groups_path and os.path.exists(groups_path):
            try:
                with open(groups_path, encoding="utf-8", errors="ignore") as fp:
                    lines = [ln.strip() for ln in fp if ln.strip()]
                # Build conv/control by taking the first entry from each line group
                conv = []
                ctrl = []
                for i, ln in enumerate(lines):
                    first_sp = next((sp.strip() for sp in ln.split(',') if sp.strip()), None)
                    if not first_sp:
                        continue
                    if i % 2 == 0:
                        conv.append(first_sp)
                    else:
                        ctrl.append(first_sp)
                if getattr(config, 'preferred_groups_combo', None) is None:
                    try:
                        config.preferred_groups_combo = (list(conv), list(ctrl))
                        config.preferred_response_matrix = ""
                    except Exception:
                        pass
            except Exception:
                conv = []
                ctrl = []

    conv = sorted(dict.fromkeys(conv))
    ctrl = sorted(dict.fromkeys(ctrl))

    # ─── Phenotype data ───────────────────────────────────────────────
    # If we parsed phenotype values from the response matrix, use them by default
    species_pheno_map: dict[str, float] | None = (resp_values if resp_values else None)
    pheno_name_map: dict[float, str] | None = None
    pheno_path = getattr(config, "species_phenotypes_file", "")
    if pheno_path and os.path.exists(pheno_path):
        try:
            species_pheno_map = {}
            with open(pheno_path) as fp:
                for line in fp:
                    parts = [p.strip() for p in line.strip().split(',') if p.strip()]
                    if len(parts) >= 2:
                        sp, val = parts[0], parts[1]
                        try:
                            species_pheno_map[sp] = float(val)
                        except ValueError:
                            continue
        except Exception:
            # Keep matrix-derived values if present; otherwise None
            species_pheno_map = species_pheno_map if species_pheno_map else None

    # Phenotype display names from config
    if hasattr(config, "pheno_name1") and hasattr(config, "pheno_name2"):
        pheno_name_map = {1.0: str(config.pheno_name1), -1.0: str(config.pheno_name2)}

    all_sites_info = None
    pss_map = {}
    if sites_path and os.path.exists(sites_path):
        try:
            df = pd.read_csv(sites_path)
            df = df[df["gene_name"] == gene]
            for _, row in df.iterrows():
                pos = int(row["position"]) - 1
                pss_map[pos] = float(row.get("pss", row.get("score", 0)))
        except Exception:
            pass

    def _save_results(self) -> None:
        # Build default filename: <align_base>__<groups_base>_site_counter_results.csv
        align_dir = getattr(self.config, 'alignments_dir', '') or ''
        align_base = os.path.basename(os.path.normpath(align_dir)) if align_dir else 'alignments'
        groups_path = getattr(self.config, 'species_groups_file', '') or ''
        groups_base = os.path.splitext(os.path.basename(groups_path))[0] if groups_path else 'groups'
        default_name = f"{align_base}__{groups_base}_site_counter_results.csv"
        # Prefer output_dir for initial location if available
        initial_dir = getattr(self.config, 'output_dir', '') or os.getcwd()
        default_path = os.path.join(initial_dir, default_name)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Site Counter Results", default_path, "CSV Files (*.csv)"
        )
        if not path:
            return
        try:
            # Prepare a copy for saving: drop redundant top_fraction columns and reorder
            df = self.results_df.copy()
            for c in ["top_fraction", "top_fraction_by_diff"]:
                if c in df.columns:
                    df = df.drop(columns=[c])
            # Desired front columns order
            front = ["gene", "num_combos_top_frac", "num_combos_top_frac_by_diff"]
            present_front = [c for c in front if c in df.columns]
            # Then metrics
            metrics = ["avg_true", "avg_control", "diff", "cs_sites_ge_4", "variable_sites", "k_pairs"]
            present_metrics = [c for c in metrics if c in df.columns]
            # Append the remaining columns (e.g., per_combo_*), preserving existing order
            remaining = [c for c in df.columns if c not in set(present_front + present_metrics)]
            ordered_cols = present_front + present_metrics + remaining
            df = df[ordered_cols]
            df.to_csv(path, index=False)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not save results:\n{e}")

    # Launch as independent top-level window (no parent) so it opens separately
    viewer = SiteViewer(
        records,
        conv,
        ctrl,
        outgroup_species or [],
        gene,
        all_sites_info,
        False,
        pss_map,
        parent=None,
        species_pheno_map=species_pheno_map,
        pheno_name_map=pheno_name_map,
    )
    # If requested (e.g., when launched from Site Counter results), default to CCS-only filter
    try:
        if prefer_ccs_filter and getattr(viewer, 'has_all_three', False):
            # Ensure filter states reflect available data, then select CCS (idx 2)
            if hasattr(viewer, 'updateFilterComboStates'):
                viewer.updateFilterComboStates()
            if hasattr(viewer, 'filter_combo') and viewer.filter_combo.model().item(2).isEnabled():
                viewer.filter_combo.setCurrentIndex(2)  # triggers rebuild via signal
    except Exception:
        # Non-fatal: if anything goes wrong, continue with default behavior
        pass
    # Ensure the viewer shows on top
    try:
        viewer.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    except Exception:
        pass
    # Show the viewer once; avoid immediate raise/activate to prevent flicker
    viewer.show()
    _open_dialogs.append(viewer)
 
def _select_combo_from_groups(parent, groups_path, current=None):
    """Open a dialog to choose a species combination from a groups file."""
    try:
        with open(groups_path, encoding="utf-8") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
    except Exception as e:
        QMessageBox.critical(parent, "Default Combo", f"Failed to read species groups file:\n{e}")
        return None
    if not lines or len(lines) % 2 != 0:
        QMessageBox.critical(parent, "Default Combo", "Species groups file has an invalid format.")
        return None
    groups = [[sp.strip() for sp in ln.split(',') if sp.strip()] for ln in lines]
    n_pairs = len(groups) // 2

    dlg = QDialog(parent)
    dlg.setWindowTitle("Select Species Combination")
    vbox = QVBoxLayout(dlg)

    conv_widgets = []
    ctrl_widgets = []
    current_conv = []
    current_ctrl = []
    if current and isinstance(current, (tuple, list)) and len(current) == 2:
        current_conv = list(current[0])
        current_ctrl = list(current[1])

    for i in range(n_pairs):
        hbox = QHBoxLayout()
        conv_label = QLabel(f"Convergent {i + 1}:")
        hbox.addWidget(conv_label)
        conv_list = QListWidget()
        conv_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        for sp in groups[2 * i]:
            conv_list.addItem(sp)
        idx = 0
        if current_conv and i < len(current_conv) and current_conv[i] in groups[2 * i]:
            idx = groups[2 * i].index(current_conv[i])
        conv_list.setCurrentRow(idx)
        if len(groups[2 * i]) == 1:
            conv_list.setEnabled(False)
        hbox.addWidget(conv_list)
        conv_widgets.append(conv_list)

        ctrl_label = QLabel(f"Control {i + 1}:")
        hbox.addWidget(ctrl_label)
        ctrl_list = QListWidget()
        ctrl_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        for sp in groups[2 * i + 1]:
            ctrl_list.addItem(sp)
        idx = 0
        if current_ctrl and i < len(current_ctrl) and current_ctrl[i] in groups[2 * i + 1]:
            idx = groups[2 * i + 1].index(current_ctrl[i])
        ctrl_list.setCurrentRow(idx)
        if len(groups[2 * i + 1]) == 1:
            ctrl_list.setEnabled(False)
        hbox.addWidget(ctrl_list)
        ctrl_widgets.append(ctrl_list)

        vbox.addLayout(hbox)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dlg)
    vbox.addWidget(buttons)
    buttons.accepted.connect(dlg.accept)
    buttons.rejected.connect(dlg.reject)

    if dlg.exec() == QDialog.Accepted:
        conv = [w.currentItem().text() for w in conv_widgets]
        ctrl = [w.currentItem().text() for w in ctrl_widgets]
        return conv, ctrl
    return None


def _get_alignment_length(gene_name: str, alignments_dir: str):
    """Return the length of the alignment for the given gene."""
    if not alignments_dir:
        return None
    path = os.path.join(alignments_dir, f"{gene_name}.fas")
    try:
        with open(path) as handle:
            for line in handle:
                if not line.startswith(">"):
                    return len(line.strip())
    except OSError:
        return None
    return None


class SpsPlotDialog(QDialog):
    """Dialog to display the SVG SPS plot."""
    def __init__(self, svg_path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SPS Prediction Plot")
        layout = QVBoxLayout(self)
        svg_widget = QSvgWidget(svg_path)
        svg_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(svg_widget)
        self.resize(680, 600)  

    @staticmethod
    def show_dialog(svg_path, parent=None):
        """Create, show, and store a reference to the dialog."""
        dialog = SpsPlotDialog(svg_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class ContinuousPlotDialog(QDialog):
    """Dialog to display the continuous phenotype plot."""

    def __init__(self, svg_path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phenotype vs SPS Plot")
        layout = QVBoxLayout(self)
        svg_widget = QSvgWidget(svg_path)
        svg_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(svg_widget)
        self.resize(680, 600)

    @staticmethod
    def show_dialog(svg_path, parent=None):
        dialog = ContinuousPlotDialog(svg_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class PredictionMetricsDialog(QDialog):
    """Display prediction metrics computed from the predictions CSV.

    Supports both binary SPS classification metrics and continuous-response
    metrics when least-squares predictions are available."""

    def __init__(self, csv_path: str, config, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Prediction Metrics")
        self._csv_path = csv_path
        self._config = config
        self._ready = False
        self._mode: str = "binary"
        self._thresholds: tuple[float, float] | None = None
        self._pearson_all: float = float("nan")
        self._last_species_stats = pd.DataFrame()

        layout = QVBoxLayout(self)
        info = QLabel("Metrics are computed from: " + os.path.basename(csv_path))
        info.setWordWrap(True)
        layout.addWidget(info)

        prog = QProgressDialog("Computing metrics…", "Cancel", 0, 40, self)
        prog.setWindowTitle("Computing Metrics")
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(200)
        step = 0
        prog.setValue(step)
        QGuiApplication.processEvents()

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to read predictions file:\n{csv_path}\n\n{e}")
            self.close()
            return
        step = self._advance_progress(prog, step)

        required = {"species", "SPS"}
        if not required.issubset(set(df.columns)):
            QMessageBox.information(
                self,
                "Unavailable",
                "Predictions CSV is missing required columns (species, SPS).",
            )
            self.close()
            return

        df_all = df.copy()
        use_continuous_predictions = bool(
            getattr(config, "use_continuous_phenotypes", False)
            or getattr(config, "response_matrices_are_continuous", False)
        )
        df_all["SPS"] = pd.to_numeric(df_all["SPS"], errors="coerce")
        if not use_continuous_predictions:
            df_all["SPS"] = df_all["SPS"].clip(lower=-1, upper=1)
        step = self._advance_progress(prog, step)

        if "input_RMSE" in df_all.columns:
            df_all["input_RMSE"] = pd.to_numeric(df_all["input_RMSE"], errors="coerce")
            df_all["RMSE_Rank"] = df_all["input_RMSE"].rank(pct=True)
        else:
            df_all["RMSE_Rank"] = float("nan")
        step = self._advance_progress(prog, step)

        if "true_phenotype" not in df.columns:
            QMessageBox.information(
                self,
                "Unavailable",
                "Prediction Metrics require a species phenotype file (true_phenotype column).",
            )
            self.close()
            return

        tp = pd.to_numeric(df_all["true_phenotype"], errors="coerce")
        df_all["tp_original"] = tp

        valid_tp = tp.dropna()
        if valid_tp.empty:
            QMessageBox.information(
                self,
                "Unavailable",
                "Prediction Metrics require at least one species with a phenotype value.",
            )
            self.close()
            return

        def _is_binary_like(values: np.ndarray) -> bool:
            allowed = (-1.0, 0.0, 1.0)
            if values.size == 0:
                return True
            for val in values:
                if not any(math.isclose(val, a, rel_tol=0.0, abs_tol=1e-6) for a in allowed):
                    return False
            return True

        is_binary_like = _is_binary_like(valid_tp.to_numpy())
        # Remember whether the true phenotypes are binary-like ({-1,0,1}).
        # We will only display Pearson correlations when phenotypes are continuous.
        self._tp_is_binary = is_binary_like
        threshold_mode = False
        if not use_continuous_predictions and not is_binary_like:
            threshold_mode = True
            dlg_thresh = PhenoThresholdDialog(valid_tp.values, parent=self)
            try:
                last_thresh = getattr(config, "_last_metrics_thresholds", None)
                if last_thresh and len(last_thresh) == 2:
                    dlg_thresh.lower_spin.setValue(float(last_thresh[0]))
                    dlg_thresh.upper_spin.setValue(float(last_thresh[1]))
            except Exception:
                pass
            if dlg_thresh.exec() != QDialog.Accepted:
                self.close()
                return
            lower = float(dlg_thresh.lower_threshold)
            upper = float(dlg_thresh.upper_threshold)
            self._thresholds = (lower, upper)
            try:
                setattr(config, "_last_metrics_thresholds", self._thresholds)
            except Exception:
                pass

            def _categorize(val: float) -> float:
                if not math.isfinite(val):
                    return float("nan")
                if val <= lower:
                    return -1.0
                if val >= upper:
                    return 1.0
                return 0.0

            df_all["tp_binary"] = df_all["tp_original"].apply(_categorize)
        else:
            def _clean_binary(val: float) -> float:
                if not math.isfinite(val):
                    return float("nan")
                for allowed in (-1.0, 0.0, 1.0):
                    if math.isclose(val, allowed, rel_tol=0.0, abs_tol=1e-6):
                        return float(allowed)
                return float("nan")

            df_all["tp_binary"] = df_all["tp_original"].apply(_clean_binary)

        if use_continuous_predictions:
            self._mode = "continuous"
        elif threshold_mode:
            self._mode = "threshold"
        else:
            self._mode = "binary"

        corr_mask = df_all["SPS"].notna() & df_all["tp_original"].notna()
        if corr_mask.sum() >= 2:
            try:
                self._pearson_all = float(
                    df_all.loc[corr_mask, "SPS"].corr(df_all.loc[corr_mask, "tp_original"])
                )
            except Exception:
                self._pearson_all = float("nan")

        def mask_all(_df: pd.DataFrame) -> pd.Series:
            return pd.Series(True, index=_df.index)

        def mask_mfs5(_df: pd.DataFrame) -> pd.Series:
            return (_df["RMSE_Rank"] < 0.05) & _df["RMSE_Rank"].notna()

        def mask_mfs10(_df: pd.DataFrame) -> pd.Series:
            return (_df["RMSE_Rank"] < 0.10) & _df["RMSE_Rank"].notna()

        def mask_mfs25(_df: pd.DataFrame) -> pd.Series:
            return (_df["RMSE_Rank"] < 0.25) & _df["RMSE_Rank"].notna()

        subsets = [
            ("All models", mask_all(df_all)),
            ("MFS bottom 5%", mask_mfs5(df_all)),
            ("MFS bottom 10%", mask_mfs10(df_all)),
            ("MFS bottom 25%", mask_mfs25(df_all)),
        ]

        extra_lines: list[str] = []
        if self._mode == "threshold" and self._thresholds is not None:
            lower, upper = self._thresholds
            extra_lines.append(
                f"Applied thresholds: ≤ {lower:.3f} classified as Control (-1), ≥ {upper:.3f} as Convergent (+1)."
            )
        # Only show Pearson correlation summary when phenotypes are continuous
        if (not getattr(self, "_tp_is_binary", False)) and math.isfinite(self._pearson_all) and not math.isnan(self._pearson_all):
            extra_lines.append(
                f"Pearson correlation between SPS and phenotype values: {self._pearson_all:.3f}."
            )
        if extra_lines:
            extra_lbl = QLabel("<br>".join(extra_lines))
            extra_lbl.setWordWrap(True)
            layout.addWidget(extra_lbl)

        if self._mode == "continuous":
            step, success = self._build_continuous_metrics(layout, df_all, subsets, prog, step)
        else:
            step, success = self._build_binary_metrics(layout, df_all, subsets, prog, step)
        if not success:
            return

        step = self._advance_progress(prog, step)
        self.resize(900, 720)
        self._ready = True
        try:
            prog.close()
        except Exception:
            pass

    def _advance_progress(self, prog: QProgressDialog, step: int, increment: int = 1) -> int:
        step += increment
        try:
            prog.setValue(min(step, prog.maximum()))
            QGuiApplication.processEvents()
        except Exception:
            pass
        return step

    @staticmethod
    def _pearson_for_mask(df_all: pd.DataFrame, mask: pd.Series) -> float:
        try:
            combined = mask & df_all["SPS"].notna() & df_all["tp_original"].notna()
        except Exception:
            return float("nan")
        work = df_all.loc[combined, ["SPS", "tp_original"]]
        if len(work) < 2:
            return float("nan")
        try:
            return float(work["SPS"].corr(work["tp_original"]))
        except Exception:
            return float("nan")

    @staticmethod
    def _compute_r2(true_vals: pd.Series | np.ndarray, pred_vals: pd.Series | np.ndarray) -> float:
        y_true = np.asarray(true_vals, dtype=float)
        y_pred = np.asarray(pred_vals, dtype=float)
        mask = np.isfinite(y_true) & np.isfinite(y_pred)
        y_true = y_true[mask]
        y_pred = y_pred[mask]
        if y_true.size < 2:
            return float("nan")
        ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
        if ss_tot <= 0.0:
            return float("nan")
        ss_res = float(np.sum((y_true - y_pred) ** 2))
        return float(1.0 - ss_res / ss_tot)

    def _format_table_value(self, column: str, value) -> str:
        if value is None:
            return ""
        if isinstance(value, (float, np.floating)):
            if np.isnan(value):
                return ""
            if column in {"Accuracy", "TPR", "TNR", "Balanced Acc", "Sign Match Frac"}:
                return f"{value:.0%}"
            if column in {"AUROC", "Pearson r", "Spearman ρ", "R²"}:
                return f"{value:.5f}"
            if column in {"RMSE", "MAE", "Bias", "Mean SPS", "Predicted Mean", "IQR"}:
                return f"{value:.5f}"
            if column == "True Phenotype":
                if self._mode == "binary":
                    return str(int(round(value)))
                return f"{value:.5f}"
            if column == "Threshold Label":
                return str(int(round(value)))
            return f"{value}"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        return str(value)

    @staticmethod
    def _acc_tpr_tnr_bal(_work: pd.DataFrame):
        labels = _work["true_phenotype"].values
        sps = _work["SPS"].values
        if len(sps) == 0:
            return 0.0, 0.0, 0.0, 0.0
        correct = (sps > 0) == (labels > 0)
        total = len(correct)
        acc = float(correct.sum()) / total if total > 0 else 0.0
        pos_mask = labels == 1
        neg_mask = labels == -1
        pos_total = int(pos_mask.sum())
        neg_total = int(neg_mask.sum())
        tpr = float((correct & pos_mask).sum()) / pos_total if pos_total > 0 else 0.0
        tnr = float((correct & neg_mask).sum()) / neg_total if neg_total > 0 else 0.0
        bal = (tpr + tnr) / 2.0 if (pos_total > 0 or neg_total > 0) else 0.0
        return acc, tpr, tnr, bal

    @staticmethod
    def _roc_auc(labels_pm1: pd.Series | np.ndarray, scores: pd.Series | np.ndarray):
        labels = np.asarray(labels_pm1)
        scores = np.asarray(scores)
        mask = np.isfinite(scores) & np.isfinite(labels)
        labels = labels[mask]
        scores = scores[mask]
        if labels.size == 0:
            return float("nan")
        labels01 = (labels > 0).astype(int)
        P = int(labels01.sum())
        N = int((1 - labels01).sum())
        if P == 0 or N == 0:
            return float("nan")
        ranks = pd.Series(scores).rank(method="average").to_numpy()
        sum_pos = ranks[labels01 == 1].sum()
        auc = (sum_pos - P * (P + 1) / 2.0) / (P * N)
        return float(auc)

    def _build_binary_metrics(
        self,
        layout: QVBoxLayout,
        df_all: pd.DataFrame,
        subsets: list[tuple[str, pd.Series]],
        prog: QProgressDialog,
        step: int,
    ) -> tuple[int, bool]:
        tp_binary = df_all["tp_binary"]
        bin_mask = tp_binary.isin([-1.0, 1.0]) & df_all["SPS"].notna()
        if bin_mask.sum() == 0:
            QMessageBox.information(
                self,
                "Unavailable",
                "No species remain after applying the binary thresholds; adjust the thresholds and try again."
                if self._mode == "threshold"
                else "Prediction Metrics are only available when labeled species (±1) are present.",
            )
            self.close()
            return step, False
        df_bin = df_all.loc[bin_mask].copy()
        df_bin["true_phenotype"] = tp_binary.loc[bin_mask].astype(int)

        # Include Pearson r column for binary/threshold metrics only if underlying phenotypes are continuous
        include_pearson = not getattr(self, "_tp_is_binary", False)

        row_cols = ["Subset", "N rows", "Accuracy", "TPR", "TNR", "Balanced Acc", "AUROC"]
        if include_pearson:
            row_cols.append("Pearson r")
        rows_summary: list[dict[str, object]] = []
        for label, mask in subsets:
            if prog.wasCanceled():
                self.close()
                return step, False
            m_on_bin = mask.reindex(df_bin.index, fill_value=False)
            work = df_bin.loc[m_on_bin]
            acc, tpr, tnr, bal = self._acc_tpr_tnr_bal(work)
            auc_rows = self._roc_auc(work["true_phenotype"], work["SPS"]) if len(work) else float("nan")
            row = {
                "Subset": label,
                "N rows": int(len(work)),
                "Accuracy": acc,
                "TPR": tpr,
                "TNR": tnr,
                "Balanced Acc": bal,
                "AUROC": auc_rows,
            }
            if include_pearson:
                row["Pearson r"] = self._pearson_for_mask(df_all, mask)
            rows_summary.append(row)
            step = self._advance_progress(prog, step)

        self.rows_table = QTableWidget(len(rows_summary), len(row_cols))
        self.rows_table.setHorizontalHeaderLabels(row_cols)
        self.rows_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.rows_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for r_idx, row in enumerate(rows_summary):
            for c_idx, col in enumerate(row_cols):
                self.rows_table.setItem(r_idx, c_idx, QTableWidgetItem(self._format_table_value(col, row.get(col))))
        self.rows_table.resizeColumnsToContents()
        layout.addWidget(QLabel("Row-level metrics"))
        layout.addWidget(self.rows_table)
        try:
            self.rows_table.resizeRowsToContents()
            vh = self.rows_table.verticalHeader()
            baseline = vh.defaultSectionSize() + 6
            for r in range(self.rows_table.rowCount()):
                self.rows_table.setRowHeight(r, max(self.rows_table.rowHeight(r), baseline))
            fix_h = self.rows_table.horizontalHeader().height() + vh.length() + 6
            self.rows_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.rows_table.setMinimumHeight(fix_h)
            self.rows_table.setMaximumHeight(fix_h)
        except Exception:
            pass

        sp_cols = ["Subset", "N species", "Accuracy", "TPR", "TNR", "Balanced Acc", "AUROC"]
        if include_pearson:
            sp_cols.append("Pearson r")
        sp_summary: list[dict[str, object]] = []
        for label, mask in subsets:
            if prog.wasCanceled():
                self.close()
                return step, False
            m_on_bin = mask.reindex(df_bin.index, fill_value=False)
            work = df_bin.loc[m_on_bin]
            sp_group = work.groupby("species", dropna=True)
            mean_sps = sp_group["SPS"].mean()
            sp_label = sp_group["true_phenotype"].first()
            n_species = int(len(mean_sps))
            if n_species > 0:
                labels_arr = sp_label.values.astype(int)
                scores_arr = mean_sps.values.astype(float)
                correct = (scores_arr > 0) == (labels_arr > 0)
                total = len(correct)
                acc = float(correct.sum()) / total if total > 0 else 0.0
                pos_mask = labels_arr == 1
                neg_mask = labels_arr == -1
                pos_total = int(pos_mask.sum())
                neg_total = int(neg_mask.sum())
                tpr = float((correct & pos_mask).sum()) / pos_total if pos_total > 0 else 0.0
                tnr = float((correct & neg_mask).sum()) / neg_total if neg_total > 0 else 0.0
                bal = (tpr + tnr) / 2.0 if (pos_total > 0 or neg_total > 0) else 0.0
                auc_species = self._roc_auc(labels_arr, scores_arr)
            else:
                acc = tpr = tnr = bal = float("nan")
                auc_species = float("nan")
            row = {
                "Subset": label,
                "N species": n_species,
                "Accuracy": acc,
                "TPR": tpr,
                "TNR": tnr,
                "Balanced Acc": bal,
                "AUROC": auc_species,
            }
            if include_pearson:
                mask_all = mask & df_all["SPS"].notna() & df_all["tp_original"].notna()
                work_all = df_all.loc[mask_all]
                pearson_species = float("nan")
                if not work_all.empty:
                    species_means = work_all.groupby("species", dropna=True).agg(
                        pred=("SPS", "mean"), true=("tp_original", "mean")
                    )
                    if len(species_means) >= 2:
                        try:
                            pearson_species = float(species_means["pred"].corr(species_means["true"]))
                        except Exception:
                            pearson_species = float("nan")
                row["Pearson r"] = pearson_species
            sp_summary.append(row)
            step = self._advance_progress(prog, step)

        self.species_summary_table = QTableWidget(len(sp_summary), len(sp_cols))
        self.species_summary_table.setHorizontalHeaderLabels(sp_cols)
        self.species_summary_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.species_summary_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for r_idx, row in enumerate(sp_summary):
            for c_idx, col in enumerate(sp_cols):
                self.species_summary_table.setItem(
                    r_idx, c_idx, QTableWidgetItem(self._format_table_value(col, row.get(col)))
                )
        self.species_summary_table.resizeColumnsToContents()
        layout.addWidget(QLabel("Species-mean metrics"))
        layout.addWidget(self.species_summary_table)
        try:
            self.species_summary_table.resizeRowsToContents()
            vh2 = self.species_summary_table.verticalHeader()
            baseline2 = vh2.defaultSectionSize() + 6
            for r in range(self.species_summary_table.rowCount()):
                self.species_summary_table.setRowHeight(
                    r, max(self.species_summary_table.rowHeight(r), baseline2)
                )
            fix_h2 = self.species_summary_table.horizontalHeader().height() + vh2.length() + 6
            self.species_summary_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.species_summary_table.setMinimumHeight(fix_h2)
            self.species_summary_table.setMaximumHeight(fix_h2)
        except Exception:
            pass

        chooser_row = QHBoxLayout()
        chooser_row.addWidget(QLabel("Per-species table for:"))
        self.subset_combo = QComboBox()
        for label, _mask in subsets:
            self.subset_combo.addItem(label)
        chooser_row.addWidget(self.subset_combo)
        chooser_row.addStretch()
        export_btn = QPushButton("Export Table…")
        chooser_row.addWidget(export_btn)
        layout.addLayout(chooser_row)

        self.species_table = QTableWidget()
        self.species_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.species_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.species_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.species_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.species_table)

        subset_map = {label: mask for (label, mask) in subsets}

        def build_species_table(which: str) -> None:
            mask = subset_map.get(which, subsets[0][1])
            work = df_all.loc[mask]
            if work.empty:
                self.species_table.clear()
                self.species_table.setRowCount(0)
                self.species_table.setColumnCount(0)
                self._last_species_stats = pd.DataFrame()
                return
            grp_all = work.groupby("species", dropna=True)
            agg = grp_all["SPS"].agg(
                mean_sps="mean",
                q25=lambda s: float(np.nanpercentile(s, 25)) if len(s) else float("nan"),
                q75=lambda s: float(np.nanpercentile(s, 75)) if len(s) else float("nan"),
                n_rows="count",
            ).reset_index()
            agg["IQR"] = agg["q75"] - agg["q25"]
            agg["true_phenotype"] = grp_all["tp_original"].first().reindex(agg["species"]).values
            agg["threshold_label"] = grp_all["tp_binary"].first().reindex(agg["species"]).values
            mask_on_bin_all = mask.reindex(df_bin.index, fill_value=False)
            work_bin_subset = df_bin.loc[mask_on_bin_all]
            grp_bin = work_bin_subset.groupby("species", dropna=True)
            try:
                match_frac = grp_bin.apply(
                    lambda g: float(((g["SPS"] > 0) == (g["true_phenotype"] > 0)).mean()) if len(g) else float("nan")
                )
            except Exception:
                match_frac = pd.Series(dtype=float)
            if not match_frac.empty:
                agg = agg.merge(match_frac.rename("sign_match_frac"), on="species", how="left")
            else:
                agg["sign_match_frac"] = float("nan")

            if self._mode == "threshold":
                headers = [
                    "Species",
                    "True Phenotype",
                    "Threshold Label",
                    "Mean SPS",
                    "IQR",
                    "N rows",
                    "Sign Match Frac",
                ]
                display_cols = [
                    "species",
                    "true_phenotype",
                    "threshold_label",
                    "mean_sps",
                    "IQR",
                    "n_rows",
                    "sign_match_frac",
                ]
            else:
                headers = ["Species", "True Phenotype", "Mean SPS", "IQR", "N rows", "Sign Match Frac"]
                display_cols = ["species", "true_phenotype", "mean_sps", "IQR", "n_rows", "sign_match_frac"]

            self.species_table.setRowCount(len(agg))
            self.species_table.setColumnCount(len(headers))
            self.species_table.setHorizontalHeaderLabels(headers)
            for r_idx, row in agg.iterrows():
                for c_idx, col in enumerate(display_cols):
                    col_name = headers[c_idx]
                    val = row.get(col)
                    self.species_table.setItem(
                        r_idx,
                        c_idx,
                        QTableWidgetItem(self._format_table_value(col_name, val)),
                    )
            self.species_table.resizeColumnsToContents()

            export_df = agg[
                ["species", "true_phenotype", "threshold_label", "mean_sps", "q25", "q75", "IQR", "n_rows", "sign_match_frac"]
            ].copy()
            if self._mode != "threshold":
                export_df["threshold_label"] = export_df["true_phenotype"]
            export_df.insert(0, "subset", which)
            export_df = export_df.rename(
                columns={
                    "species": "species",
                    "true_phenotype": "true_phenotype",
                    "threshold_label": "threshold_label",
                    "mean_sps": "mean_sps",
                    "q25": "q25",
                    "q75": "q75",
                    "IQR": "iqr",
                    "n_rows": "n_rows",
                    "sign_match_frac": "sign_match_frac",
                }
            )
            self._last_species_stats = export_df

        def on_subset_changed(_idx: int) -> None:
            build_species_table(self.subset_combo.currentText())

        def export_current_table() -> None:
            preds_dir = os.path.dirname(self._csv_path)
            preds_base = os.path.splitext(os.path.basename(self._csv_path))[0]
            subset_label = self.subset_combo.currentText()
            token_map = {
                "All models": "all_models",
                "MFS bottom 5%": "5pct_mfs",
                "MFS bottom 10%": "10pct_mfs",
                "MFS bottom 25%": "25pct_mfs",
            }
            token = token_map.get(subset_label, "all_models")
            default_name = f"per-species_predictions_{token}_{preds_base}.csv"
            default_path = os.path.join(preds_dir, default_name)
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Species Metrics", default_path, "CSV Files (*.csv)"
            )
            if not path:
                return
            try:
                stats = getattr(self, "_last_species_stats", None)
                if stats is None or stats.empty:
                    return
                stats.to_csv(path, index=False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not save CSV:\n{e}")

        self.subset_combo.currentIndexChanged.connect(on_subset_changed)
        export_btn.clicked.connect(export_current_table)
        build_species_table(self.subset_combo.currentText())
        step = self._advance_progress(prog, step)
        return step, True

    def _build_continuous_metrics(
        self,
        layout: QVBoxLayout,
        df_all: pd.DataFrame,
        subsets: list[tuple[str, pd.Series]],
        prog: QProgressDialog,
        step: int,
    ) -> tuple[int, bool]:
        valid_mask = df_all["SPS"].notna() & df_all["tp_original"].notna()
        if valid_mask.sum() == 0:
            QMessageBox.information(
                self,
                "Unavailable",
                "No predictions contain both SPS values and continuous phenotypes.",
            )
            self.close()
            return step, False

        row_cols = ["Subset", "N rows", "Pearson r", "Spearman ρ", "RMSE", "MAE", "R²"]
        rows_summary: list[dict[str, object]] = []
        for label, mask in subsets:
            if prog.wasCanceled():
                self.close()
                return step, False
            mask_rows = mask & df_all["SPS"].notna() & df_all["tp_original"].notna()
            work = df_all.loc[mask_rows]
            n_rows = int(len(work))
            if n_rows >= 2:
                pearson = float(work["SPS"].corr(work["tp_original"]))
                spearman = float(work["SPS"].corr(work["tp_original"], method="spearman"))
            else:
                pearson = spearman = float("nan")
            diff = work["SPS"] - work["tp_original"]
            mae = float(np.abs(diff).mean()) if n_rows else float("nan")
            rmse = float(np.sqrt(np.mean(np.square(diff)))) if n_rows else float("nan")
            r2 = self._compute_r2(work["tp_original"], work["SPS"])
            rows_summary.append(
                {
                    "Subset": label,
                    "N rows": n_rows,
                    "Pearson r": pearson,
                    "Spearman ρ": spearman,
                    "RMSE": rmse,
                    "MAE": mae,
                    "R²": r2,
                }
            )
            step = self._advance_progress(prog, step)

        self.rows_table = QTableWidget(len(rows_summary), len(row_cols))
        self.rows_table.setHorizontalHeaderLabels(row_cols)
        self.rows_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.rows_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for r_idx, row in enumerate(rows_summary):
            for c_idx, col in enumerate(row_cols):
                self.rows_table.setItem(
                    r_idx, c_idx, QTableWidgetItem(self._format_table_value(col, row.get(col)))
                )
        self.rows_table.resizeColumnsToContents()
        layout.addWidget(QLabel("Row-level continuous metrics"))
        layout.addWidget(self.rows_table)
        try:
            self.rows_table.resizeRowsToContents()
            vh = self.rows_table.verticalHeader()
            baseline = vh.defaultSectionSize() + 6
            for r in range(self.rows_table.rowCount()):
                self.rows_table.setRowHeight(r, max(self.rows_table.rowHeight(r), baseline))
            fix_h = self.rows_table.horizontalHeader().height() + vh.length() + 6
            self.rows_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.rows_table.setMinimumHeight(fix_h)
            self.rows_table.setMaximumHeight(fix_h)
        except Exception:
            pass

        sp_cols = ["Subset", "N species", "Pearson r", "Spearman ρ", "RMSE", "MAE", "R²"]
        sp_summary: list[dict[str, object]] = []
        for label, mask in subsets:
            if prog.wasCanceled():
                self.close()
                return step, False
            mask_rows = mask & df_all["SPS"].notna() & df_all["tp_original"].notna()
            work = df_all.loc[mask_rows]
            grp = work.groupby("species", dropna=True)
            means = grp.agg(pred=("SPS", "mean"), true=("tp_original", "mean"))
            n_species = int(len(means))
            if n_species >= 2:
                pearson = float(means["pred"].corr(means["true"]))
                spearman = float(means["pred"].corr(means["true"], method="spearman"))
            else:
                pearson = spearman = float("nan")
            err = means["pred"] - means["true"]
            mae = float(np.abs(err).mean()) if n_species else float("nan")
            rmse = float(np.sqrt(np.mean(np.square(err)))) if n_species else float("nan")
            r2 = self._compute_r2(means["true"], means["pred"])
            sp_summary.append(
                {
                    "Subset": label,
                    "N species": n_species,
                    "Pearson r": pearson,
                    "Spearman ρ": spearman,
                    "RMSE": rmse,
                    "MAE": mae,
                    "R²": r2,
                }
            )
            step = self._advance_progress(prog, step)

        self.species_summary_table = QTableWidget(len(sp_summary), len(sp_cols))
        self.species_summary_table.setHorizontalHeaderLabels(sp_cols)
        self.species_summary_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.species_summary_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for r_idx, row in enumerate(sp_summary):
            for c_idx, col in enumerate(sp_cols):
                self.species_summary_table.setItem(
                    r_idx, c_idx, QTableWidgetItem(self._format_table_value(col, row.get(col)))
                )
        self.species_summary_table.resizeColumnsToContents()
        layout.addWidget(QLabel("Species-mean continuous metrics"))
        layout.addWidget(self.species_summary_table)
        try:
            self.species_summary_table.resizeRowsToContents()
            vh2 = self.species_summary_table.verticalHeader()
            baseline2 = vh2.defaultSectionSize() + 6
            for r in range(self.species_summary_table.rowCount()):
                self.species_summary_table.setRowHeight(
                    r, max(self.species_summary_table.rowHeight(r), baseline2)
                )
            fix_h2 = self.species_summary_table.horizontalHeader().height() + vh2.length() + 6
            self.species_summary_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.species_summary_table.setMinimumHeight(fix_h2)
            self.species_summary_table.setMaximumHeight(fix_h2)
        except Exception:
            pass

        chooser_row = QHBoxLayout()
        chooser_row.addWidget(QLabel("Per-species table for:"))
        self.subset_combo = QComboBox()
        for label, _mask in subsets:
            self.subset_combo.addItem(label)
        chooser_row.addWidget(self.subset_combo)
        chooser_row.addStretch()
        export_btn = QPushButton("Export Table…")
        chooser_row.addWidget(export_btn)
        layout.addLayout(chooser_row)

        self.species_table = QTableWidget()
        self.species_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.species_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.species_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.species_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.species_table)

        subset_map = {label: mask for (label, mask) in subsets}

        def build_species_table(which: str) -> None:
            mask = subset_map.get(which, subsets[0][1])
            mask_rows = mask & df_all["SPS"].notna() & df_all["tp_original"].notna()
            work = df_all.loc[mask_rows]
            if work.empty:
                self.species_table.clear()
                self.species_table.setRowCount(0)
                self.species_table.setColumnCount(0)
                self._last_species_stats = pd.DataFrame()
                return
            grp = work.groupby("species", dropna=True)
            stats = grp["SPS"].agg(
                predicted_mean="mean",
                q25=lambda s: float(np.nanpercentile(s, 25)) if len(s) else float("nan"),
                q75=lambda s: float(np.nanpercentile(s, 75)) if len(s) else float("nan"),
                n_rows="count",
            ).reset_index()
            stats["IQR"] = stats["q75"] - stats["q25"]
            stats["true_phenotype"] = grp["tp_original"].mean().reindex(stats["species"]).values
            stats["bias"] = stats["predicted_mean"] - stats["true_phenotype"]

            def _mae_series(g: pd.DataFrame) -> float:
                diff = g["SPS"] - g["tp_original"]
                return float(np.abs(diff).mean()) if len(diff) else float("nan")

            def _rmse_series(g: pd.DataFrame) -> float:
                diff = g["SPS"] - g["tp_original"]
                return float(np.sqrt(np.mean(np.square(diff)))) if len(diff) else float("nan")

            mae = grp.apply(_mae_series)
            rmse = grp.apply(_rmse_series)
            stats = stats.merge(mae.rename("mae"), left_on="species", right_index=True, how="left")
            stats = stats.merge(rmse.rename("rmse"), left_on="species", right_index=True, how="left")

            headers = [
                "Species",
                "True Phenotype",
                "Predicted Mean",
                "Bias",
                "MAE",
                "RMSE",
                "IQR",
                "N rows",
            ]
            display_cols = [
                "species",
                "true_phenotype",
                "predicted_mean",
                "bias",
                "mae",
                "rmse",
                "IQR",
                "n_rows",
            ]
            self.species_table.setRowCount(len(stats))
            self.species_table.setColumnCount(len(headers))
            self.species_table.setHorizontalHeaderLabels(headers)
            for r_idx, row in stats.iterrows():
                for c_idx, col in enumerate(display_cols):
                    col_name = headers[c_idx]
                    val = row.get(col)
                    self.species_table.setItem(
                        r_idx,
                        c_idx,
                        QTableWidgetItem(self._format_table_value(col_name, val)),
                    )
            self.species_table.resizeColumnsToContents()

            export_df = stats[
                ["species", "true_phenotype", "predicted_mean", "bias", "mae", "rmse", "q25", "q75", "IQR", "n_rows"]
            ].copy()
            export_df.insert(0, "subset", which)
            export_df = export_df.rename(
                columns={
                    "species": "species",
                    "true_phenotype": "true_phenotype",
                    "predicted_mean": "predicted_mean",
                    "bias": "bias",
                    "mae": "mae",
                    "rmse": "rmse",
                    "q25": "q25",
                    "q75": "q75",
                    "IQR": "iqr",
                    "n_rows": "n_rows",
                }
            )
            self._last_species_stats = export_df

        def on_subset_changed(_idx: int) -> None:
            build_species_table(self.subset_combo.currentText())

        def export_current_table() -> None:
            preds_dir = os.path.dirname(self._csv_path)
            preds_base = os.path.splitext(os.path.basename(self._csv_path))[0]
            subset_label = self.subset_combo.currentText()
            token_map = {
                "All models": "all_models",
                "MFS bottom 5%": "5pct_mfs",
                "MFS bottom 10%": "10pct_mfs",
                "MFS bottom 25%": "25pct_mfs",
            }
            token = token_map.get(subset_label, "all_models")
            default_name = f"per-species_predictions_{token}_{preds_base}.csv"
            default_path = os.path.join(preds_dir, default_name)
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Species Metrics", default_path, "CSV Files (*.csv)"
            )
            if not path:
                return
            try:
                stats = getattr(self, "_last_species_stats", None)
                if stats is None or stats.empty:
                    return
                stats.to_csv(path, index=False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not save CSV:\n{e}")

        self.subset_combo.currentIndexChanged.connect(on_subset_changed)
        export_btn.clicked.connect(export_current_table)
        build_species_table(self.subset_combo.currentText())
        step = self._advance_progress(prog, step)
        return step, True

    @staticmethod
    def show_dialog(csv_path: str, config, parent=None):
        dialog = PredictionMetricsDialog(csv_path, config, parent)
        if dialog._ready:
            dialog.show()
            _open_dialogs.append(dialog)

class GeneRanksDialog(QWidget):
    """Dialog to display top gene ranks and optional selected sites."""

    def __init__(self, dataframe, config, sites_path=None, parent=None):
        super().__init__(parent)
        # Make this a true normal top-level window (not Dialog) for stacking
        try:
            base_flags = (
                Qt.WindowType.Window
                | Qt.WindowType.WindowTitleHint
                | Qt.WindowType.WindowSystemMenuHint
                | Qt.WindowType.WindowMinimizeButtonHint
                | Qt.WindowType.WindowMaximizeButtonHint
                | Qt.WindowType.WindowCloseButtonHint
            )
            self.setWindowFlags(base_flags)
        except Exception:
            pass
        self.setWindowTitle("Top Gene Ranks")
        layout = QVBoxLayout(self)

        self.config = config
        self.sites_path = sites_path
        self._selected_groups_combo = getattr(config, 'preferred_groups_combo', None)
        self._selected_response_matrix = getattr(config, 'preferred_response_matrix', '')

        df = dataframe

        # Clean and rename columns for readability
        rename_map = {
            'gene_name': 'Gene',
            'highest_gss': 'Highest GSS',
            'highest_ever_gss': 'Highest Ever GSS',
            'best_rank': 'Best Rank',
            'best_ever_rank': 'Best Ever Rank',
            'num_selected_sites': 'Number of Selected Sites',
            'num_combos_ranked': 'Number of Combinations Ranked',
            'num_combos_ranked_top': 'Number Ranked in Top',
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        # Round and format columns
        if 'Highest GSS' in df.columns:
            df['Highest GSS'] = df['Highest GSS'].round(5)
        if 'Highest Ever GSS' in df.columns:
            df['Highest Ever GSS'] = pd.to_numeric(df['Highest Ever GSS'], errors='coerce').round(6)
        # Ensure rank columns are numeric while safely handling placeholders like 'None'
        for col in ['Best Rank', 'Best Ever Rank']:
            if col in df.columns:
                # Convert values to numeric, setting invalid parsing (e.g. 'None', '') to NaN
                numeric_vals = pd.to_numeric(df[col], errors='coerce')
                # Replace NaN with empty string to keep the table display tidy
                df[col] = numeric_vals.apply(lambda x: '' if pd.isna(x) else int(x))

        self.has_sites = bool(sites_path and os.path.exists(sites_path))

        msg = "Double-click a row to open the Site Viewer."
        if self.has_sites:
            msg += " Use the expand/collapse arrows to reveal selected sites for each gene."
        help_label = QLabel(msg)
        # Keep the help text on a single line and let it expand horizontally
        help_label.setWordWrap(False)
        help_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        # Place help text and a 'Set Combo' button on the same row
        header_row = QHBoxLayout()
        header_row.addWidget(help_label)
        header_row.setStretchFactor(help_label, 1)
        header_row.addStretch()
        set_combo_btn = QPushButton("Set Combo")
        set_combo_btn.setToolTip("Choose which response matrix combo to use when opening Site Viewer windows.")
        set_combo_btn.clicked.connect(self._open_combo_picker)
        header_row.addWidget(set_combo_btn)
        layout.addLayout(header_row)

        if self.has_sites:
            self._init_tree_view(layout, df)
        else:
            self._init_table_view(layout, df)

        # Set a wider initial size to better accommodate the content
        self.resize(1200, 800)
        
        # Ensure the window is scrollable if content is too wide
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        
        # Set a minimum size to prevent the window from being too small
        self.setMinimumSize(1000, 600)

        # Enable copying selected rows with Ctrl/Cmd+C
        try:
            sc = QShortcut(QKeySequence.Copy, self)
            sc.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            sc.activated.connect(self._copy_selected_rows_to_clipboard)
        except Exception:
            pass

    def _open_site_viewer(self, row_or_item, _column: int) -> None:
        """Launch :class:`SiteViewer` for the selected gene."""
        # When using the tree view the signal provides a ``QTreeWidgetItem``
        # whereas the table view provides the row index.  Handle both cases.
        if self.has_sites:
            # If called from ``QTreeWidget.itemDoubleClicked`` we are passed the
            # item directly.  The older implementation expected a row index,
            # which caused a crash because ``topLevelItem`` requires an int.  To
            # remain backwards compatible we accept either form.
            if isinstance(row_or_item, QTreeWidgetItem):
                item = row_or_item
            else:
                item = self.tree.topLevelItem(int(row_or_item)) if row_or_item is not None else None
            # Only act on top-level gene rows, not the child site rows.
            if item is None or item.parent() is not None:
                return
            gene = item.text(0)
        else:
            row = int(row_or_item)
            gene_item = self.table.item(row, self.gene_col)
            if gene_item is None:
                return
            gene = gene_item.text()
        try:
            _launch_site_viewer(gene, self.config, self.sites_path, parent=self)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Site Viewer:\n{e}")

    def _open_combo_picker(self):
        """Open a dialog allowing the user to pick species groups or a response matrix."""
        groups_path = getattr(self.config, 'species_groups_file', '') or ''
        if groups_path and os.path.exists(groups_path):
            current = getattr(self, '_selected_groups_combo', None)
            res = _select_combo_from_groups(self, groups_path, current)
            if res:
                try:
                    setattr(self, '_selected_groups_combo', res)
                    self.config.preferred_groups_combo = res
                    self.config.preferred_response_matrix = ""
                except Exception:
                    pass
            return

        response_dir = getattr(self.config, 'response_dir', '') or ''
        if not response_dir:
            base = os.path.basename(getattr(self.config, 'species_groups_file', '')).replace('.txt', '')
            if base and getattr(self.config, 'output_dir', ''):
                response_dir = os.path.join(self.config.output_dir, f"{base}_response_matrices")
        if not response_dir or not os.path.isdir(response_dir):
            QMessageBox.information(self, "No Response Matrices", "No response matrix directory was found.")
            return

        files = sorted([f for f in os.listdir(response_dir) if f.endswith('.txt')])
        if not files:
            QMessageBox.information(self, "No Response Matrices", "No .txt response matrices found in the directory.")
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Select Response Combo")
        vbox = QVBoxLayout(dlg)
        combo = QComboBox(dlg)
        combo.addItems(files)
        try:
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
            fm = QFontMetrics(combo.font())
            longest = max(files, key=len)
            text_w = fm.horizontalAdvance(longest)
            min_w = min(900, text_w + 120)
            combo.setMinimumWidth(min_w)
            dlg.setMinimumWidth(min_w + 60)
        except Exception:
            dlg.resize(520, dlg.sizeHint().height())
        current = getattr(self, '_selected_response_matrix', '') or ''
        if current and current in files:
            combo.setCurrentIndex(files.index(current))
        vbox.addWidget(combo)

        conv_label = QLabel("Convergent species:")
        conv_list = QLabel("")
        conv_list.setWordWrap(True)
        ctrl_label = QLabel("Control species:")
        ctrl_list = QLabel("")
        ctrl_list.setWordWrap(True)
        vbox.addWidget(conv_label)
        vbox.addWidget(conv_list)
        vbox.addWidget(ctrl_label)
        vbox.addWidget(ctrl_list)

        def parse_species_lists(fname: str):
            path = os.path.join(response_dir, fname)
            conv, ctrl = [], []
            try:
                with open(path, encoding='utf-8', errors='ignore') as f:
                    for idx, raw in enumerate(f):
                        line = raw.strip()
                        if not line:
                            continue
                        parts = line.split()
                        if len(parts) < 2:
                            continue
                        sp = parts[0]
                        if idx % 2 == 0:
                            conv.append(sp)
                        else:
                            ctrl.append(sp)
            except Exception:
                pass
            return conv, ctrl

        def refresh_lists():
            fname = combo.currentText()
            conv, ctrl = parse_species_lists(fname)
            conv_list.setText("\n".join(conv) if conv else "(none)")
            ctrl_list.setText("\n".join(ctrl) if ctrl else "(none)")

        combo.currentIndexChanged.connect(lambda _i: refresh_lists())
        refresh_lists()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dlg)
        vbox.addWidget(buttons)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)

        if dlg.exec() == QDialog.Accepted:
            selected = combo.currentText()
            try:
                setattr(self, '_selected_response_matrix', selected)
                self.config.preferred_response_matrix = selected
                self.config.preferred_groups_combo = None
            except Exception:
                pass

    def _init_table_view(self, layout: QVBoxLayout, df: pd.DataFrame) -> None:
        """Initialize the simple table view (no selected sites)."""
        self.table = QTableWidget(len(df.index), len(df.columns))
        self.table.setHorizontalHeaderLabels(list(df.columns))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        # Enable selecting arbitrary cells and multiple ranges for copying
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setSortingEnabled(False)
        self.gene_col = list(df.columns).index('Gene') if 'Gene' in df.columns else 0

        for r_idx, (_, row) in enumerate(df.iterrows()):
            for c_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r_idx, c_idx, item)

        self.table.resizeColumnsToContents()
        layout.addWidget(self.table)
        self.table.cellDoubleClicked.connect(self._open_site_viewer)

    def _copy_selected_rows_to_clipboard(self):
        """Copy selection to the clipboard as TSV.

        - In tree view (with sites): copy full top-level gene rows for selected items.
        - In table view: copy selected rectangular cell blocks (supports multiple ranges).
        """
        try:
            lines = []
            if getattr(self, 'has_sites', False):
                # Copy top-level gene rows from the tree
                if not hasattr(self, 'tree') or self.tree is None:
                    return
                items = self.tree.selectedItems()
                if not items:
                    return
                # Collect unique top-level items
                top_items = []
                seen = set()
                for it in items:
                    root = it
                    while root.parent() is not None:
                        root = root.parent()
                    if id(root) not in seen:
                        seen.add(id(root))
                        top_items.append(root)
                cols = self.tree.columnCount()
                for it in top_items:
                    vals = [it.text(c) for c in range(cols)]
                    lines.append('\t'.join(vals))
            else:
                # Copy selected rectangular cell ranges from the table
                if not hasattr(self, 'table') or self.table is None:
                    return
                ranges = self.table.selectedRanges()
                if not ranges:
                    return
                blocks = []
                for r in ranges:
                    block_lines = []
                    for row in range(r.topRow(), r.bottomRow() + 1):
                        vals = []
                        for col in range(r.leftColumn(), r.rightColumn() + 1):
                            item = self.table.item(row, col)
                            vals.append(item.text() if item is not None else '')
                        block_lines.append('\t'.join(vals))
                    blocks.append('\n'.join(block_lines))
                lines = ['\n'.join(blocks)]
            if lines:
                QGuiApplication.clipboard().setText('\n'.join(lines))
        except Exception:
            pass

    def _init_tree_view(self, layout: QVBoxLayout, df: pd.DataFrame) -> None:
        """Initialize the tree view with selected sites information."""
        try:
            sites_df = pd.read_csv(self.sites_path)
        except Exception:
            # Fallback to table view if parsing fails
            self.has_sites = False
            self._init_table_view(layout, df)
            return

        self.tree = QTreeWidget()

        rank_columns = [c for c in df.columns if c != 'Gene']
        headers = ['Gene', 'Length', 'Map', 'Position', 'Max PSS'] + rank_columns
        self.tree.setColumnCount(len(headers))
        self.tree.setHeaderLabels(headers)
        self.tree.setUniformRowHeights(True)
        self.tree.header().setSectionResizeMode(2, QHeaderView.ResizeMode.Fixed)
        self.tree.setColumnWidth(2, 520)

        align_dir = getattr(self.config, 'alignments_dir', '')

        # Compute a global maximum score across the displayed genes (up to 200)
        global_max_score = None
        try:
            display_genes = set(df['Gene']) if 'Gene' in df.columns else set()
            if display_genes:
                sub = sites_df[sites_df['gene_name'].isin(display_genes)]
                score_col = None
                if 'pss' in sub.columns and not sub['pss'].isnull().all():
                    score_col = 'pss'
                elif 'score' in sub.columns:
                    score_col = 'score'
                if score_col is not None and not sub.empty:
                    mx = sub[score_col].max()
                    if pd.notna(mx) and float(mx) > 0:
                        global_max_score = float(mx)
        except Exception:
            global_max_score = None

        for _, row in df.iterrows():
            gene = row['Gene']
            gene_sites = sites_df[sites_df['gene_name'] == gene]
            gene_sites = gene_sites[gene_sites['pss'] != 0] if 'pss' in gene_sites.columns else gene_sites
            positions = gene_sites['position'].tolist()
            scores = (
                gene_sites['pss'].tolist()
                if 'pss' in gene_sites.columns and not gene_sites['pss'].isnull().all()
                else gene_sites['score'].tolist() if 'score' in gene_sites.columns else []
            )
            # Pre-compute score range for consistent coloring
            min_score = min(scores) if scores else None
            max_score = max(scores) if scores else None
            if min_score is not None and max_score == min_score:
                max_score += 1e-9
            length = _get_alignment_length(gene, align_dir) or (max(positions) if positions else 1)

            gene_values = [gene, str(length), '', '', ''] + [str(row[col]) for col in rank_columns]
            parent_item = QTreeWidgetItem(gene_values)
            self.tree.addTopLevelItem(parent_item)
            self.tree.setItemWidget(parent_item, 2, ProteinMapWidget(length, positions, scores, score_scale_max=global_max_score))
            parent_item.setExpanded(False)

            for _, srow in gene_sites.iterrows():
                pos = str(srow['position'])
                pss_val = srow['pss'] if 'pss' in srow else srow.get('score', 0)
                if pss_val == 0:
                    continue
                pss = f"{pss_val:.5f}"
                child_vals = ['', '', '', pos, pss] + ['' for _ in rank_columns]
                child = QTreeWidgetItem(child_vals)
                # Apply the same color scaling used by ProteinMapWidget
                if min_score is not None and max_score is not None:
                    rel = (pss_val - min_score) / (max_score - min_score)
                    brightness = 0.2 + 0.8 * rel
                else:
                    brightness = 1.0
                base_color = QColor("#4CAF50")
                r = int(base_color.red() * brightness)
                g = int(base_color.green() * brightness)
                b = int(base_color.blue() * brightness)
                brush = QBrush(QColor(r, g, b))
                child.setForeground(4, brush)
                parent_item.addChild(child)

        self.tree.resizeColumnToContents(0)
        self.tree.resizeColumnToContents(1)
        self.tree.resizeColumnToContents(3)
        self.tree.resizeColumnToContents(4)

        layout.addWidget(self.tree)
        self.tree.itemDoubleClicked.connect(self._open_site_viewer)

    @staticmethod
    def show_dialog(csv_path, config, sites_path=None, parent=None):
        """Load data, then create, show, and store a reference to the dialog."""
        try:
            # Read file as-is and keep the existing row order (already sorted by the CLI).
            df = pd.read_csv(csv_path).head(200)
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to read gene ranks file:\n{csv_path}\n\n{e}")
            return

        dialog = GeneRanksDialog(df, config, sites_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class SelectedSitesDialog(QDialog):
    """(Deprecated) Dialog to display selected sites grouped by gene."""
    def __init__(self, dataframe, gene_order, alignments_dir, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Selected Sites")
        layout = QVBoxLayout(self)
        
        df = dataframe

        tree = QTreeWidget()
        tree.setColumnCount(5)
        tree.setHeaderLabels(["Gene", "Length", "Map", "Position", "Position Sparsity Score"])
        tree.setUniformRowHeights(True)

        # Compute a global maximum score across the displayed genes
        global_max_score = None
        try:
            sub = df[df["gene_name"].isin(gene_order)] if "gene_name" in df.columns else df
            score_col = None
            if "pss" in sub.columns and not sub["pss"].isnull().all():
                score_col = "pss"
            elif "score" in sub.columns:
                score_col = "score"
            if score_col is not None and not sub.empty:
                mx = sub[score_col].max()
                if pd.notna(mx) and float(mx) > 0:
                    global_max_score = float(mx)
        except Exception:
            global_max_score = None

        for gene in gene_order:
            gene_data = df[df["gene_name"] == gene]
            # Filter out rows with zero PSS
            gene_data = gene_data[gene_data["pss"] != 0] if "pss" in gene_data.columns else gene_data
            if gene_data.empty:
                continue

            positions = gene_data["position"].tolist()
            scores = (
                gene_data["pss"].tolist()
                if "pss" in gene_data.columns and not gene_data["pss"].isnull().all()
                else gene_data["score"].tolist() if "score" in gene_data.columns else []
            )
            length = _get_alignment_length(gene, alignments_dir) or (max(positions) if positions else 1)

            parent_item = QTreeWidgetItem([gene, str(length), "", "", ""])
            tree.addTopLevelItem(parent_item)
            tree.setItemWidget(parent_item, 2, ProteinMapWidget(length, positions, scores, score_scale_max=global_max_score))
            parent_item.setExpanded(False)

            for _, row in gene_data.iterrows():
                pos = str(row["position"])
                pss_val = row["pss"] if "pss" in row else row.get("score", 0)
                if pss_val == 0:
                    continue  # skip zero PSS
                pss = f"{pss_val:.5f}"
                child = QTreeWidgetItem(["", "", "", pos, pss])
                parent_item.addChild(child)
        
        tree.resizeColumnToContents(0)
        tree.resizeColumnToContents(1)
        tree.resizeColumnToContents(3)
        tree.resizeColumnToContents(4)

        layout.addWidget(tree)
        self.resize(800, 600)

    @staticmethod
    def show_dialog(csv_path, gene_ranks_path, alignments_dir, parent=None):
        """Load data, then create, show, and store a reference to the dialog."""
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to read selected sites file:\n{csv_path}\n\n{e}")
            return

        # Determine gene display order
        gene_order = []
        if gene_ranks_path and os.path.exists(gene_ranks_path):
            try:
                ranks_df = pd.read_csv(gene_ranks_path)
                gene_order = list(ranks_df["gene_name"].head(200))
            except Exception:
                pass  # Fallback to order in sites file

        if not gene_order:
            gene_order = list(dict.fromkeys(df["gene_name"]))[:200]

        dialog = SelectedSitesDialog(df, gene_order, alignments_dir, parent)
        dialog.show()
        _open_dialogs.append(dialog)

class SiteCounterResultsDialog(QWidget):
    """Display site counter gene rankings."""

    def __init__(self, results, config, outgroup, parent=None):
        super().__init__(parent)
        # Make this a true normal top-level window (not Dialog) for stacking
        try:
            base_flags = (
                Qt.WindowType.Window
                | Qt.WindowType.WindowTitleHint
                | Qt.WindowType.WindowSystemMenuHint
                | Qt.WindowType.WindowMinimizeButtonHint
                | Qt.WindowType.WindowMaximizeButtonHint
                | Qt.WindowType.WindowCloseButtonHint
            )
            self.setWindowFlags(base_flags)
        except Exception:
            pass
        self.setWindowTitle("Site Counter Results")
        self.config = config
        self.outgroup = outgroup
        self._selected_groups_combo = getattr(config, 'preferred_groups_combo', None)
        self._selected_response_matrix = getattr(config, 'preferred_response_matrix', '')
        layout = QVBoxLayout(self)

        header = QHBoxLayout()
        label = QLabel("Double-click a gene to open the Site Viewer.")
        header.addWidget(label)
        header.addStretch()
        # Buttons: Save Results first (focused), then Set Combo to its right
        self.save_btn = QPushButton("Save Results")
        self.save_btn.clicked.connect(self._save_results)
        header.addWidget(self.save_btn)
        self.set_combo_btn = QPushButton("Set Combo")
        self.set_combo_btn.setToolTip("Choose a response matrix (if provided) or a species-groups combo to use when opening the Site Viewer from these results.")
        self.set_combo_btn.clicked.connect(self._open_combo_picker)
        header.addWidget(self.set_combo_btn)
        layout.addLayout(header)

        self.table = QTableWidget()
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        # Allow selecting arbitrary cells (not just rows) and multiple ranges
        self.table.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.doubleClicked.connect(lambda idx: self._open_site_viewer(idx.row()))
        # Enable clickable header sorting; NumericItem ensures numeric sort for numeric columns
        self.table.setSortingEnabled(True)
        try:
            hh = self.table.horizontalHeader()
            hh.setSortIndicatorShown(True)
            hh.setSectionsClickable(True)
        except Exception:
            pass
        layout.addWidget(self.table)

        self.results_df = pd.DataFrame(results)

        # Configure columns dynamically: include combo-based columns only when present
        has_combo_rank_true = 'num_combos_top_frac' in self.results_df.columns
        has_combo_rank_ratio = 'num_combos_top_frac_by_ratio' in self.results_df.columns
        has_combo_rank_diff = 'num_combos_top_frac_by_diff' in self.results_df.columns
        # Single-combo mode: no combo-ranking columns present. Compute a per-gene ratio
        # consistent with multi-combo logic: ratio = true / (control + 1), where
        # control is derived from per_combo_diff when available (ctrl = max(0, t - d)).
        single_combo_mode = not (has_combo_rank_true or has_combo_rank_ratio or has_combo_rank_diff)
        if single_combo_mode:
            def _compute_single_ratio(row):
                try:
                    tvals = row.get('per_combo_true', None)
                    dvals = row.get('per_combo_diff', None)
                    t = None
                    d = None
                    if isinstance(tvals, (list, tuple)) and len(tvals) >= 1:
                        t = tvals[0]
                    if isinstance(dvals, (list, tuple)) and len(dvals) >= 1:
                        d = dvals[0]
                    if t is None:
                        at = row.get('avg_true', None)
                        if at is None or pd.isna(at):
                            return float('nan')
                        t = float(at)
                    else:
                        t = float(t)
                    ctrl = 0.0
                    if d is not None:
                        try:
                            d = float(d)
                            ctrl = max(0.0, t - d)
                        except Exception:
                            ctrl = 0.0
                    return t / (ctrl + 1.0)
                except Exception:
                    return float('nan')
            try:
                self.results_df['ratio'] = self.results_df.apply(_compute_single_ratio, axis=1)
            except Exception:
                pass
        sort_col_idx = None
        if has_combo_rank_true or has_combo_rank_ratio or has_combo_rank_diff:
            # Build headers with 'by Ratio' then 'by Diff', then plain True after those
            headers = ["Gene"]
            idx_ratio_hdr = None
            idx_diff_hdr = None
            idx_true_hdr = None
            if has_combo_rank_ratio:
                idx_ratio_hdr = len(headers)
                headers.append("Combos in Top % by Ratio")
            if has_combo_rank_diff:
                idx_diff_hdr = len(headers)
                headers.append("Combos in Top % by Diff")
            if has_combo_rank_true:
                idx_true_hdr = len(headers)
                headers.append("Combos in Top %")
            # Core metrics
            headers += [
                "Avg True Convergence",
                "Avg Control Convergence",
                "Avg True - Control",
                "CS ≥ 4 Sites",
                "Variable Sites",
            ]
            # Inject percentages into headers if available and record default sort column
            try:
                pct = None
                if has_combo_rank_true:
                    tf_series = self.results_df.get('top_fraction')
                    if tf_series is not None and not tf_series.dropna().empty:
                        pct = int(round(float(tf_series.dropna().iloc[0]) * 100))
                        if idx_true_hdr is not None:
                            headers[idx_true_hdr] = f"Combos in Top {pct}%"
                # Ratio uses its own stored fraction but falls back to the True's if not present
                if has_combo_rank_ratio:
                    tfr = self.results_df.get('top_fraction_by_ratio')
                    pct_r = None
                    if tfr is not None and not tfr.dropna().empty:
                        pct_r = int(round(float(tfr.dropna().iloc[0]) * 100))
                    elif pct is not None:
                        pct_r = pct
                    if pct_r is not None and idx_ratio_hdr is not None:
                        headers[idx_ratio_hdr] = f"Combos in Top {pct_r}% by Ratio"
                if has_combo_rank_diff:
                    tf_series_d = self.results_df.get('top_fraction_by_diff')
                    pct_d = None
                    if tf_series_d is not None and not tf_series_d.dropna().empty:
                        pct_d = int(round(float(tf_series_d.dropna().iloc[0]) * 100))
                    elif pct is not None:
                        pct_d = pct
                    if pct_d is not None and idx_diff_hdr is not None:
                        headers[idx_diff_hdr] = f"Combos in Top {pct_d}% by Diff"
                # Default sort column -> prefer Ratio if present (more robust), else True, else Diff
                if idx_ratio_hdr is not None:
                    sort_col_idx = idx_ratio_hdr
                elif idx_diff_hdr is not None:
                    sort_col_idx = idx_diff_hdr
                elif idx_true_hdr is not None:
                    sort_col_idx = idx_true_hdr
            except Exception:
                pass
            self.table.setColumnCount(len(headers))
            self.table.setHorizontalHeaderLabels(headers)
        else:
            headers = [
                "Gene",
                "True/(Control+1) Ratio",
                "Avg True Convergence",
                "Avg Control Convergence",
                "Avg True - Control",
                "CS ≥ 4 Sites",
                "Variable Sites",
            ]
            self.table.setColumnCount(7)
            self.table.setHorizontalHeaderLabels(headers)
            # Fallback sort on Ratio
            sort_col_idx = 1
        # Show all rows where average true convergence > 0 (include all, no 200 cap)
        if 'avg_true' in self.results_df.columns:
            filtered = self.results_df[self.results_df['avg_true'] > 0]
            top = filtered if not filtered.empty else self.results_df
        else:
            top = self.results_df
        self.table.setRowCount(len(top))
        def _fmt_num(v) -> str:
            try:
                x = float(v)
            except Exception:
                return str(v)
            if math.isclose(x, round(x), rel_tol=0.0, abs_tol=1e-9):
                return str(int(round(x)))
            # Up to 3 decimals, strip trailing zeros and trailing dot
            s = f"{x:.3f}".rstrip('0').rstrip('.')
            return s

        for row_idx, (_, row) in enumerate(top.iterrows()):
            # Gene column – regular text item
            self.table.setItem(row_idx, 0, QTableWidgetItem(str(row["gene"])) )
            col = 1
            # Combos by Ratio first (if present)
            if has_combo_rank_ratio:
                combos_top_r = row.get('num_combos_top_frac_by_ratio', None)
                if combos_top_r is None or (isinstance(combos_top_r, float) and pd.isna(combos_top_r)):
                    display_ct_r = ''
                    combos_top_val_r = float('nan')
                else:
                    try:
                        combos_top_val_r = float(combos_top_r)
                    except Exception:
                        combos_top_val_r = float('nan')
                    display_ct_r = _fmt_num(combos_top_val_r) if combos_top_val_r == combos_top_val_r else ''
                self.table.setItem(row_idx, col, NumericItem(combos_top_val_r, display_ct_r))
                col += 1
            # Combos by Diff next (if present)
            if has_combo_rank_diff:
                combos_top_d = row.get('num_combos_top_frac_by_diff', None)
                if combos_top_d is None or (isinstance(combos_top_d, float) and pd.isna(combos_top_d)):
                    display_ct_d = ''
                    combos_top_val_d = float('nan')
                else:
                    try:
                        combos_top_val_d = float(combos_top_d)
                    except Exception:
                        combos_top_val_d = float('nan')
                    display_ct_d = _fmt_num(combos_top_val_d) if combos_top_val_d == combos_top_val_d else ''
                self.table.setItem(row_idx, col, NumericItem(combos_top_val_d, display_ct_d))
                col += 1
            # Plain "Combos in Top %" (by true) last, if present
            if has_combo_rank_true:
                combos_top = row.get('num_combos_top_frac', None)
                if combos_top is None or (isinstance(combos_top, float) and pd.isna(combos_top)):
                    display_ct = ''
                    combos_top_val = float('nan')
                else:
                    try:
                        combos_top_val = float(combos_top)
                    except Exception:
                        combos_top_val = float('nan')
                    display_ct = _fmt_num(combos_top_val) if combos_top_val == combos_top_val else ''
                self.table.setItem(row_idx, col, NumericItem(combos_top_val, display_ct))
                col += 1
            # Insert single-combo ratio first if applicable
            if single_combo_mode:
                rv = row.get('ratio', None)
                if rv is None or (isinstance(rv, float) and pd.isna(rv)):
                    display_ratio = ''
                    rvf = float('nan')
                else:
                    try:
                        rvf = float(rv)
                    except Exception:
                        rvf = float('nan')
                    display_ratio = _fmt_num(rvf) if rvf == rvf else ''
                self.table.setItem(row_idx, col, NumericItem(rvf, display_ratio))
                col += 1
            # Core metrics
            v1 = float(row['avg_true']) if pd.notna(row['avg_true']) else float('nan')
            v2 = float(row['avg_control']) if pd.notna(row['avg_control']) else float('nan')
            v3 = float(row['diff']) if pd.notna(row['diff']) else float('nan')
            self.table.setItem(row_idx, col,   NumericItem(v1, _fmt_num(v1))); col += 1
            self.table.setItem(row_idx, col,   NumericItem(v2, _fmt_num(v2))); col += 1
            self.table.setItem(row_idx, col,   NumericItem(v3, _fmt_num(v3))); col += 1
            # CS ≥ 4 sites
            cs_sites = float(row.get('cs_sites_ge_4', float('nan')))
            self.table.setItem(row_idx, col,   NumericItem(cs_sites, _fmt_num(cs_sites))); col += 1
            # Variable Sites at far right
            v_sites = float(row['variable_sites']) if pd.notna(row['variable_sites']) else float('nan')
            self.table.setItem(row_idx, col,   NumericItem(v_sites, _fmt_num(v_sites)))
        self.table.resizeColumnsToContents()
        self.resize(800, 600)
        # Make Save Results the default and focused button
        try:
            self.save_btn.setAutoDefault(True)
            self.save_btn.setDefault(True)
            self.save_btn.setFocus()
        except Exception:
            pass

        # Apply default sort: by Diff-based combos count if present, else by Avg True
        try:
            if sort_col_idx is not None:
                from PySide6.QtCore import Qt as _Qt
                self.table.sortItems(sort_col_idx, _Qt.SortOrder.DescendingOrder)
        except Exception:
            pass

        # Enable copying selected rows with Ctrl/Cmd+C
        try:
            sc = QShortcut(QKeySequence.Copy, self)
            sc.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            sc.activated.connect(self._copy_selected_rows_to_clipboard)
        except Exception:
            pass

    def _copy_selected_rows_to_clipboard(self):
        """Copy selected rectangular cell ranges from the site counter table as TSV."""
        try:
            if not hasattr(self, 'table') or self.table is None:
                return
            ranges = self.table.selectedRanges()
            if not ranges:
                return
            blocks = []
            for r in ranges:
                block_lines = []
                for row in range(r.topRow(), r.bottomRow() + 1):
                    vals = []
                    for col in range(r.leftColumn(), r.rightColumn() + 1):
                        item = self.table.item(row, col)
                        vals.append(item.text() if item is not None else '')
                    block_lines.append('\t'.join(vals))
                blocks.append('\n'.join(block_lines))
            if blocks:
                QGuiApplication.clipboard().setText('\n'.join(blocks))
        except Exception:
            pass

    def _open_combo_picker(self) -> None:
        """Open a dialog to pick the default combo for Site Viewer."""
        groups_path = getattr(self.config, 'species_groups_file', '') or ''
        if groups_path and os.path.exists(groups_path):
            current = getattr(self, '_selected_groups_combo', None)
            res = _select_combo_from_groups(self, groups_path, current)
            if res:
                try:
                    setattr(self, '_selected_groups_combo', res)
                    self.config.preferred_groups_combo = res
                    self.config.preferred_response_matrix = ""
                except Exception:
                    pass
            return

        response_dir = getattr(self.config, 'response_dir', '') or ''
        if response_dir and os.path.isdir(response_dir):
            files = sorted([f for f in os.listdir(response_dir) if f.endswith('.txt')])
            if not files:
                QMessageBox.information(self, "Default Combo", "No response matrices found in the selected directory.")
                return
            dlg = QDialog(self)
            dlg.setWindowTitle("Select Response Combo")
            vbox = QVBoxLayout(dlg)
            combo = QComboBox(dlg)
            combo.addItems(files)
            try:
                combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
                fm = QFontMetrics(combo.font())
                longest = max(files, key=len)
                text_w = fm.horizontalAdvance(longest)
                combo.setMinimumWidth(min(max(text_w + 40, 240), 640))
            except Exception:
                pass

            vbox.addWidget(QLabel("Select Combo:"))
            vbox.addWidget(combo)

            def parse_species_lists(fname: str):
                path = os.path.join(response_dir, fname)
                conv, ctrl = [], []
                try:
                    with open(path, encoding='utf-8', errors='ignore') as f:
                        for idx, raw in enumerate(f):
                            line = raw.strip()
                            if not line:
                                continue
                            parts = line.split()
                            if len(parts) < 2:
                                continue
                            sp = parts[0]
                            if idx % 2 == 0:
                                conv.append(sp)
                            else:
                                ctrl.append(sp)
                except Exception:
                    pass
                return conv, ctrl

            conv_list = QLabel()
            ctrl_list = QLabel()
            conv_list.setTextInteractionFlags(Qt.TextSelectableByMouse)
            ctrl_list.setTextInteractionFlags(Qt.TextSelectableByMouse)
            vbox.addWidget(QLabel("Convergent (trait-positive) species:"))
            vbox.addWidget(conv_list)
            vbox.addWidget(QLabel("Control (trait-negative) species:"))
            vbox.addWidget(ctrl_list)

            def refresh_lists():
                fname = combo.currentText()
                conv, ctrl = parse_species_lists(fname)
                conv_list.setText("\n".join(conv) if conv else "(none)")
                ctrl_list.setText("\n".join(ctrl) if ctrl else "(none)")

            combo.currentIndexChanged.connect(lambda _i: refresh_lists())
            refresh_lists()

            buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dlg)
            vbox.addWidget(buttons)
            buttons.accepted.connect(dlg.accept)
            buttons.rejected.connect(dlg.reject)

            if dlg.exec() == QDialog.Accepted:
                selected = combo.currentText()
                try:
                    setattr(self, '_selected_response_matrix', selected)
                    self.config.preferred_response_matrix = selected
                    self.config.preferred_groups_combo = None
                except Exception:
                    pass
            return

        QMessageBox.information(self, "Default Combo", "No species groups file or response directory available.")

    def _open_site_viewer(self, row: int) -> None:
        gene_item = self.table.item(row, 0)
        if gene_item is None:
            return
        gene = gene_item.text()
        try:
            _launch_site_viewer(
                gene,
                self.config,
                None,
                parent=self,
                outgroup_species=[self.outgroup],
                prefer_ccs_filter=True,
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Site Viewer:\n{e}")

    def _save_results(self) -> None:
        # Build default filename: <align_base>__<groups_base>_site_counter_results.csv
        align_dir = getattr(self.config, 'alignments_dir', '') or ''
        align_base = os.path.basename(os.path.normpath(align_dir)) if align_dir else 'alignments'
        groups_path = getattr(self.config, 'species_groups_file', '') or ''
        groups_base = os.path.splitext(os.path.basename(groups_path))[0] if groups_path else 'groups'
        default_name = f"{align_base}__{groups_base}_site_counter_results.csv"
        # Prefer output_dir for initial location if available
        initial_dir = getattr(self.config, 'output_dir', '') or os.getcwd()
        default_path = os.path.join(initial_dir, default_name)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Site Counter Results", default_path, "CSV Files (*.csv)"
        )
        if not path:
            return
        try:
            # Drop extremely wide per-combo arrays from the export to keep CSV manageable
            df = self.results_df.copy()
            for col in ("per_combo_true", "per_combo_diff"):
                if col in df.columns:
                    df = df.drop(columns=[col])
            # Reorder columns: Gene, then combo-rank columns (if present), then core metrics, then the rest
            present_front = [c for c in ["gene", "ratio", "num_combos_top_frac_by_ratio", "num_combos_top_frac_by_diff", "num_combos_top_frac"] if c in df.columns]
            metrics = [
                "avg_true",
                "avg_control",
                "diff",
                "cs_sites_ge_4",
                "variable_sites",
                "k_pairs",
            ]
            present_metrics = [c for c in metrics if c in df.columns]
            remaining = [c for c in df.columns if c not in set(present_front + present_metrics)]
            ordered_cols = present_front + present_metrics + remaining
            df = df[ordered_cols]
            df.to_csv(path, index=False)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not save results:\n{e}")

    @staticmethod
    def show_results(results, config, outgroup, parent=None):
        dialog = SiteCounterResultsDialog(results, config, outgroup, parent)
        dialog.show()
        _open_dialogs.append(dialog)


# Backward compatibility alias for older imports.
FastScanResultsDialog = SiteCounterResultsDialog
