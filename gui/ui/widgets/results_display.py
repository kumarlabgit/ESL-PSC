"""
Dialogs for displaying ESL-PSC analysis results.
"""
from __future__ import annotations
import os

import math
# Use shared FASTA reader
from gui.core.fasta_io import read_fasta
from esl_psc_cli.deletion_canceler import (
    parse_species_groups as cli_parse_species_groups,
)
import pandas as pd
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QBrush, QFontMetrics, QKeySequence, QGuiApplication, QShortcut
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QSizePolicy, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QTreeWidget, QTreeWidgetItem, QMessageBox,
    QLabel, QHeaderView, QPushButton, QHBoxLayout,
    QComboBox, QDialogButtonBox, QFileDialog
)
from PySide6.QtSvgWidgets import QSvgWidget

from gui.ui.widgets.protein_map import ProteinMapWidget
from gui.ui.widgets.site_viewer import SiteViewer

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
    # This enables Fast Scan to open SiteViewer with a user-chosen species grouping
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
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Fast Scan Results", "fast_scan_results.csv", "CSV Files (*.csv)"
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
    # Ensure the viewer shows on top
    try:
        viewer.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    except Exception:
        pass
    # Raise the results dialog first so it stays above the wizard, then show
    # the viewer and bring it to the front.
    try:
        if parent is not None:
            # Ensure native window exists then raise without taking focus
            _ = parent.winId()
            parent.raise_()
    except Exception:
        pass
    viewer.show()
    try:
        viewer.raise_()
        viewer.activateWindow()
    except Exception:
        pass

    # After the event loop processes the new window, reassert the stacking
    # by raising the viewer again. Do this twice to handle platform/WM quirks.
    def _reassert_order():
        try:
            viewer.raise_()
            viewer.activateWindow()
        except Exception:
            pass
    try:
        QTimer.singleShot(0, _reassert_order)
        QTimer.singleShot(50, _reassert_order)
    except Exception:
        pass

    # When the viewer closes, bring the launching dialog back to the front so
    # it doesn't appear to "disappear" behind the wizard window.
    if parent is not None:
        def _refocus_parent():
            try:
                parent.show()
                # On some platforms, toggling WindowActive ensures focus
                try:
                    parent.raise_()
                    parent.activateWindow()
                except Exception:
                    pass
            except Exception:
                pass
        try:
            viewer.destroyed.connect(lambda _obj=None: _refocus_parent())
        except Exception:
            pass
    _open_dialogs.append(viewer)

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


class GeneRanksDialog(QDialog):
    """Dialog to display top gene ranks and optional selected sites."""

    def __init__(self, dataframe, config, sites_path=None, parent=None):
        super().__init__(parent)
        # Make this a true top-level window (behaves like a normal app window)
        try:
            self.setWindowFlag(Qt.WindowType.Window, True)
        except Exception:
            pass
        self.setWindowTitle("Top Gene Ranks")
        layout = QVBoxLayout(self)

        self.config = config
        self.sites_path = sites_path

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
        """Open a dialog allowing the user to pick a response matrix combo and preview its species."""
        response_dir = getattr(self.config, 'response_dir', '') or ''
        if not response_dir:
            # Attempt to derive from species_groups_file fallback
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

        # Build the dialog
        dlg = QDialog(self)
        dlg.setWindowTitle("Select Response Combo")
        vbox = QVBoxLayout(dlg)
        combo = QComboBox(dlg)
        combo.addItems(files)
        # Adjust sizing so the dialog/title comfortably fits
        try:
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
            fm = QFontMetrics(combo.font())
            longest = max(files, key=len)
            text_w = fm.horizontalAdvance(longest)
            min_w = min(900, text_w + 120)  # padding for arrow/margins
            combo.setMinimumWidth(min_w)
            dlg.setMinimumWidth(min_w + 60)
        except Exception:
            # Reasonable fallback width
            dlg.resize(520, dlg.sizeHint().height())
        # Preselect current dialog selection if present
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
                        parts = line.split()  # tab or whitespace
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
            # Persist selection for subsequent Site Viewer openings on this dialog instance
            try:
                setattr(self, '_selected_response_matrix', selected)
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

class FastScanResultsDialog(QDialog):
    """Display fast scan gene rankings."""

    def __init__(self, results, config, outgroup, parent=None):
        super().__init__(parent)
        # Make this a true top-level window (behaves like a normal app window)
        try:
            self.setWindowFlag(Qt.WindowType.Window, True)
        except Exception:
            pass
        self.setWindowTitle("Fast Scan Results")
        self.config = config
        self.outgroup = outgroup
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
        has_combo_rank_diff = 'num_combos_top_frac_by_diff' in self.results_df.columns
        sort_col_idx = None
        if has_combo_rank_true or has_combo_rank_diff:
            # Build headers with 'by True' immediately after Gene, then optional 'by Diff'
            headers = ["Gene"]
            idx_true_hdr = None
            idx_diff_hdr = None
            if has_combo_rank_true:
                idx_true_hdr = len(headers)
                headers.append("Combos in Top %")
            if has_combo_rank_diff:
                idx_diff_hdr = len(headers)
                headers.append("Combos in Top % by Diff")
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
                if has_combo_rank_diff:
                    tf_series_d = self.results_df.get('top_fraction_by_diff')
                    pct_d = None
                    if tf_series_d is not None and not tf_series_d.dropna().empty:
                        pct_d = int(round(float(tf_series_d.dropna().iloc[0]) * 100))
                    elif pct is not None:
                        pct_d = pct
                    if pct_d is not None and idx_diff_hdr is not None:
                        headers[idx_diff_hdr] = f"Combos in Top {pct_d}% by Diff"
                # Default sort column -> by True if present, else by Diff if present
                if idx_true_hdr is not None:
                    sort_col_idx = idx_true_hdr
                elif idx_diff_hdr is not None:
                    sort_col_idx = idx_diff_hdr
            except Exception:
                pass
            self.table.setColumnCount(len(headers))
            self.table.setHorizontalHeaderLabels(headers)
        else:
            headers = [
                "Gene",
                "Avg True Convergence",
                "Avg Control Convergence",
                "Avg True - Control",
                "CS ≥ 4 Sites",
                "Variable Sites",
            ]
            self.table.setColumnCount(6)
            self.table.setHorizontalHeaderLabels(headers)
            # Fallback sort on Avg True
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
            # Combos by True immediately after Gene, if present
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
        """Copy selected rectangular cell ranges from the fast scan table as TSV."""
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
        """Open a dialog to pick the default combo for Site Viewer.

        Priority:
        - If a response_dir is set (and exists), allow picking a response matrix file
          (same behavior as the main ESL-PSC results display dialog).
        - Otherwise, parse combos from the species groups file and allow picking one
          of the computed Cartesian-product combos.
        """
        # Try response_dir first
        response_dir = getattr(self.config, 'response_dir', '') or ''
        if response_dir and os.path.isdir(response_dir):
            files = sorted([f for f in os.listdir(response_dir) if f.endswith('.txt')])
            if not files:
                QMessageBox.information(self, "Default Combo", "No response matrices found in the selected directory.")
                return
            from PySide6.QtWidgets import QDialogButtonBox
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

            # Place the selection dropdown at the top, matching the main results dialog UX
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
                except Exception:
                    pass
            return

        # Otherwise, parse combos from species groups file
        groups_path = getattr(self.config, 'species_groups_file', '') or ''
        if not groups_path or not os.path.exists(groups_path):
            QMessageBox.information(self, "Default Combo", "No species groups file or response directory available.")
            return
        try:
            raw_combos = cli_parse_species_groups(groups_path)
        except Exception as e:
            QMessageBox.critical(self, "Default Combo", f"Failed to parse species groups file:\n{e}")
            return
        if not raw_combos:
            QMessageBox.information(self, "Default Combo", "No combos found in species groups file.")
            return

        # Build display items and map to conv/ctrl
        display_items = []
        conv_ctrl_map = {}
        for i, tup in enumerate(raw_combos):
            picks = list(tup)
            conv = [picks[j] for j in range(0, len(picks), 2)]
            ctrl = [picks[j] for j in range(1, len(picks), 2)]
            key = f"combo_{i}"
            display_items.append(key)
            conv_ctrl_map[key] = (conv, ctrl)

        from PySide6.QtWidgets import QDialogButtonBox
        dlg = QDialog(self)
        dlg.setWindowTitle("Select Groups Combo")
        vbox = QVBoxLayout(dlg)
        combo = QComboBox(dlg)
        combo.addItems(display_items)
        try:
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
            fm = QFontMetrics(combo.font())
            longest = max(display_items, key=len)
            text_w = fm.horizontalAdvance(longest)
            combo.setMinimumWidth(min(max(text_w + 40, 240), 640))
        except Exception:
            pass

        # Place the selection dropdown at the top, matching the main results dialog UX
        vbox.addWidget(QLabel("Select Combo:"))
        vbox.addWidget(combo)

        conv_list = QLabel()
        ctrl_list = QLabel()
        conv_list.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ctrl_list.setTextInteractionFlags(Qt.TextSelectableByMouse)
        vbox.addWidget(QLabel("Convergent (trait-positive) species:"))
        vbox.addWidget(conv_list)
        vbox.addWidget(QLabel("Control (trait-negative) species:"))
        vbox.addWidget(ctrl_list)

        def refresh_lists():
            key = combo.currentText()
            conv, ctrl = conv_ctrl_map.get(key, ([], []))
            conv_list.setText("\n".join(conv) if conv else "(none)")
            ctrl_list.setText("\n".join(ctrl) if ctrl else "(none)")

        combo.currentIndexChanged.connect(lambda _i: refresh_lists())
        refresh_lists()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dlg)
        vbox.addWidget(buttons)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)

        if dlg.exec() == QDialog.Accepted:
            key = combo.currentText()
            conv, ctrl = conv_ctrl_map.get(key, ([], []))
            try:
                setattr(self, '_selected_groups_combo', (conv, ctrl))
            except Exception:
                pass

    def _open_site_viewer(self, row: int) -> None:
        gene_item = self.table.item(row, 0)
        if gene_item is None:
            return
        gene = gene_item.text()
        try:
            _launch_site_viewer(gene, self.config, None, parent=self, outgroup_species=[self.outgroup])
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Site Viewer:\n{e}")

    def _save_results(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Fast Scan Results", "fast_scan_results.csv", "CSV Files (*.csv)"
        )
        if path:
            try:
                self.results_df.to_csv(path, index=False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not save results:\n{e}")

    @staticmethod
    def show_results(results, config, outgroup, parent=None):
        dialog = FastScanResultsDialog(results, config, outgroup, parent)
        dialog.show()
        _open_dialogs.append(dialog)