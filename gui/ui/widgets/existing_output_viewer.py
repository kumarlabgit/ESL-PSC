"""Dialog for viewing results from an existing ESL-PSC run.

Allows the user to load a previously completed analysis output directory (must
contain a ``checkpoint/command.json`` file) and provides quick-access buttons to
view the generated SPS plot and gene-rank table using the same viewer dialogs
employed by the Run page.
"""
from __future__ import annotations

import os
from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QMessageBox,
    QFileDialog,
)

from gui.core.config import ESLConfig
from gui.ui.widgets.results_display import SpsPlotDialog, GeneRanksDialog

__all__ = [
    "select_and_show_existing_output",
]


class _ExistingOutputDialog(QDialog):
    """Internal dialog that displays basic run information & result buttons."""

    def __init__(
        self,
        cfg: ESLConfig,
        output_dir: str,
        parent=None,
    ) -> None:
        super().__init__(parent)

        self._cfg = cfg
        self._output_dir = output_dir

        self.setWindowTitle("Existing ESL-PSC Output Viewer")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        layout = QVBoxLayout(self)

        # Heading
        heading = QLabel(f"<h3>Viewing output for: <code>{os.path.basename(output_dir)}</code></h3>")
        heading.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(heading)

        # ─── Key run information ──────────────────────────────────────────
        info_lines = [
            f"<b>Alignments directory:</b> {cfg.alignments_dir or 'N/A'}",
            f"<b>Species groups file:</b> {cfg.species_groups_file or 'N/A'}",
            f"<b>Output directory:</b> {cfg.output_dir}",
            f"<b>Output base name:</b> {cfg.output_file_base_name}",
        ]
        lbl = QLabel("<br>".join(info_lines))
        lbl.setTextFormat(Qt.TextFormat.RichText)
        lbl.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse | Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        layout.addWidget(lbl)

        # ─── Result buttons ───────────────────────────────────────────────
        btn_layout = QHBoxLayout()
        layout.addLayout(btn_layout)
        btn_layout.addStretch()

        self._sps_path: Optional[str] = self._find_sps_plot()
        self._gene_ranks_path: Optional[str] = self._find_gene_ranks()
        self._sites_path: Optional[str] = self._find_selected_sites()

        self._sps_btn = QPushButton("Show SPS Plot")
        if self._sps_path:
            self._sps_btn.setToolTip("Open the SPS prediction plot from this run")
        else:
            self._sps_btn.setEnabled(False)
            self._sps_btn.setToolTip("No SPS plot found in this output folder")
            from PySide6.QtWidgets import QGraphicsOpacityEffect
            eff = QGraphicsOpacityEffect(self._sps_btn)
            eff.setOpacity(0.4)
            self._sps_btn.setGraphicsEffect(eff)
        self._sps_btn.clicked.connect(self._show_sps)
        btn_layout.addWidget(self._sps_btn)

        self._gene_btn = QPushButton("Show Gene Ranks")
        if self._gene_ranks_path:
            self._gene_btn.setToolTip("Open the gene-rank table from this run")
        else:
            self._gene_btn.setEnabled(False)
            self._gene_btn.setToolTip("No gene-ranks file found in this output folder")
            from PySide6.QtWidgets import QGraphicsOpacityEffect
            eff = QGraphicsOpacityEffect(self._gene_btn)
            eff.setOpacity(0.4)
            self._gene_btn.setGraphicsEffect(eff)
        self._gene_btn.clicked.connect(self._show_gene_ranks)
        btn_layout.addWidget(self._gene_btn)

        btn_layout.addStretch()

        self.resize(520, 180)

    # ------------------------------------------------------------------
    # File helpers
    # ------------------------------------------------------------------
    def _find_sps_plot(self) -> Optional[str]:
        base = self._cfg.output_file_base_name
        if not base or getattr(self._cfg, "no_pred_output", False):
            return None
        cand = os.path.join(self._output_dir, f"{base}_pred_sps_plot.svg")
        return cand if os.path.exists(cand) else None

    def _find_gene_ranks(self) -> Optional[str]:
        base = self._cfg.output_file_base_name
        if not base or getattr(self._cfg, "no_genes_output", False):
            return None
        cand = os.path.join(self._output_dir, f"{base}_gene_ranks.csv")
        return cand if os.path.exists(cand) else None

    def _find_selected_sites(self) -> Optional[str]:
        base = self._cfg.output_file_base_name
        cand = os.path.join(self._output_dir, f"{base}_selected_sites.csv")
        return cand if os.path.exists(cand) else None

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------
    def _show_sps(self):  # pragma: no cover – UI slot
        if self._sps_path:
            SpsPlotDialog.show_dialog(self._sps_path, parent=self)
        else:
            QMessageBox.information(self, "Unavailable", "No SPS plot found for this run.")

    def _show_gene_ranks(self):  # pragma: no cover – UI slot
        if self._gene_ranks_path:
            GeneRanksDialog.show_dialog(
                self._gene_ranks_path,
                self._cfg,
                self._sites_path,
                parent=self,
            )
        else:
            QMessageBox.information(self, "Unavailable", "No gene-ranks file found for this run.")


# ──────────────────────────────────────────────────────────────────────────────
# Public helper
# ──────────────────────────────────────────────────────────────────────────────

def select_and_show_existing_output(parent=None):
    """Prompt user to choose a finished ESL-PSC output directory and launch the viewer.

    A *valid* folder must contain either a gene-ranks CSV or an SPS plot **and**
    a corresponding ``*_run_config.txt`` file written at run start. The
    run-config is parsed to rebuild an :class:`~gui.core.config.ESLConfig``
    instance for displaying metadata.
    """
    dir_path = QFileDialog.getExistingDirectory(
        parent,
        "Select ESL-PSC Output Folder",
        os.getcwd(),
    )
    if not dir_path:
        return  # user cancelled

    # Locate the *first* run-config file in the folder (pattern: *_run_config.txt)
    run_cfg_files = [f for f in os.listdir(dir_path) if f.endswith("_run_config.txt")]
    if not run_cfg_files:
        QMessageBox.warning(
            parent,
            "Invalid Folder",
            "The selected directory does not contain a run-config file ( *_run_config.txt ).",
        )
        return

    run_cfg_path = os.path.join(dir_path, run_cfg_files[0])


    # ------------------------------------------------------------------
    # Parse run-config – each line mirrors CLI args ("--flag [value]")
    # ------------------------------------------------------------------
    import shlex
    arg_dict = {}
    try:
        with open(run_cfg_path, "r", encoding="utf-8") as fh:
            for raw in fh:
                raw = raw.strip()
                if not raw:
                    continue
                parts = shlex.split(raw)
                flag = parts[0]
                name = flag.lstrip("-")
                if len(parts) == 1:
                    arg_dict[name] = True  # store_true flag
                else:
                    # Join rest back with space for paths that may contain spaces
                    arg_dict[name] = " ".join(parts[1:])
    except Exception as exc:
        QMessageBox.critical(parent, "Error", f"Failed to read run-config file:\n{exc}")
        return

    # Build ESLConfig from stored command args where possible
    cfg = ESLConfig()
    for key, val in arg_dict.items():
        if hasattr(cfg, key):
            setattr(cfg, key, val)
    # Fallback defaults
    if not cfg.output_dir:
        cfg.output_dir = dir_path

    # Validate presence of at least one displayable output
    base = cfg.output_file_base_name or ""
    has_gene = os.path.exists(os.path.join(dir_path, f"{base}_gene_ranks.csv"))
    has_sps = os.path.exists(os.path.join(dir_path, f"{base}_pred_sps_plot.svg"))
    if not (has_gene or has_sps):
        QMessageBox.warning(
            parent,
            "No Viewable Output",
            "The selected directory does not contain a gene-ranks CSV or SPS plot to display.",
        )
        return

    dlg = _ExistingOutputDialog(cfg, cfg.output_dir, parent=parent)
    dlg.show()
