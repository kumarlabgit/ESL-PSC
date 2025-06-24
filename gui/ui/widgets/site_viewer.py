# widgets/site_viewer.py
from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush, QFont
from PyQt6.QtWidgets import (
    QTableWidget, QTableWidgetItem, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QSlider, QComboBox, QLabel, QPushButton,
    QCheckBox, QAbstractItemView, QMessageBox, QMenu
)

from gui.ui.widgets.histogram_canvas import HistogramCanvas
from gui.constants import ZAPPO_STATIC_COLORS



class SiteViewer(QWidget):
    """Main alignment-inspection widget with persistent outgroup memory between
    SiteViewer windows within the same GUI session."""

    # ------------------------------------------------------------------
    # Class-level cache to remember the last outgroup species list chosen by
    # the user. This is shared across *all* SiteViewer instances created in the
    # current Python process.
    # ------------------------------------------------------------------
    REMEMBERED_OUTGROUP: List[str] = []

    """
    Main alignment-inspection widget.
    """

    def __init__(
        self,
        records,                # List[tuple[id, seq]]
        convergent_species: List[str],
        control_species:    List[str],
        outgroup_species:   List[str],
        gene_name: str | None = None,
        all_sites_info:     List[Dict[str, Any]] | None = None,
        show_all_by_default: bool = False,
        pss_scores: Dict[int, float] | None = None,
        parent=None,
        # Optional phenotype information
        species_pheno_map: Dict[str, int] | None = None,
        pheno_name_map: Dict[int, str] | None = None,
    ) -> None:
        super().__init__(parent)

        # ------------------------------------------------------------------
        # If caller did not supply an outgroup list, attempt to re-apply the
        # remembered species from a previous SiteViewer window (if any).
        # ------------------------------------------------------------------
        if SiteViewer.REMEMBERED_OUTGROUP:
            # Move any remembered species found in convergent/control into the
            # outgroup list for this new viewer.
            # Try to add each remembered species to the outgroup if it exists in this
            # alignment and is not already there.
            all_species_this_alignment = [r[0] for r in records]
            for sp in SiteViewer.REMEMBERED_OUTGROUP:
                if sp not in all_species_this_alignment:
                    continue  # Species absent from this alignment

                # Remove it from any other explicit group lists first
                if sp in convergent_species:
                    convergent_species.remove(sp)
                    outgroup_species.append(sp)
                elif sp in control_species:
                    control_species.remove(sp)
                    outgroup_species.append(sp)
                else:
                    # Species was not explicitly assigned to any group yet
                    if sp not in outgroup_species:
                        outgroup_species.append(sp)
            # Ensure deterministic ordering
            outgroup_species.sort()

        # ------------------------------------------------------------------
        # store the minimal bits up front
        self.records              = records
        self.convergent_species  = sorted(convergent_species)
        self.control_species     = sorted(control_species)
        self.outgroup_species    = sorted(outgroup_species)
        self.pss_scores          = pss_scores or {}
        self.gene_name           = gene_name or ""

        # ─── Phenotype maps ───────────────────────────────────────────────
        # Mapping of species → phenotype value (1/-1)
        self.species_pheno_map: Dict[str, int] = species_pheno_map or {}
        # Mapping of phenotype value (1/-1) → display name (e.g. "C4", "C3")
        default_pheno_map = {1: "1", -1: "-1"}
        # Build phenotype display names mapping avoiding Python 3.9+ dict union
        self.pheno_name_map: Dict[int, str] = default_pheno_map.copy()
        if pheno_name_map:
            self.pheno_name_map.update(pheno_name_map)

        # Cache the phenotype value considered "convergent" (defaults to +1)
        self.convergent_pheno_value: int = 1

        if all_sites_info is None:
            all_sites_info = [
                {
                    'position': i,
                    'converge_degree': 0,
                    'is_ccs': False,
                    'conv_same_as_ctrl': False,
                }
                for i in range(len(records[0][1]) if records else 0)
            ]
        self.all_sites_info      = all_sites_info

        self.all_species = sorted([r[0] for r in records])
        self.species_ids = [r[0] for r in records]
        self.sequences = [r[1] for r in records]
        self.seq_length = len(self.sequences[0]) if self.sequences else 0

        self.scores = [s['converge_degree'] for s in self.all_sites_info]
        self.min_score = min(self.scores) if self.scores else 0
        self.max_score = max(self.scores) if self.scores else 0

        # Use constant default threshold of 2; will be clamped later if needed
        self.default_threshold = 2
        self.current_threshold = self.default_threshold
        # Flag to ensure we only auto-align the threshold once
        self._initial_threshold_aligned = False

        self.default_sort_mode = "position"
        # Always show all species; the checkbox will be hidden but kept in the layout
        self.show_all_species = True

        # Additional controls if we have all three groups
        self.only_ccs = False
        self.hide_control_convergence = False
        # Filter for ESL-PSC selected sites
        self.only_selected = False
        self.has_selected_sites = bool(self.pss_scores)

        # Are all 3 groups non-empty?
        self.has_all_three = bool(self.convergent_species) and bool(self.control_species) and bool(self.outgroup_species)

        self._syncing_horizontal_splitters = False

        self.initUI()
        self.rebuildTables()
        # Adjust vertical splitter sizes based on table content
        self._adjustVerticalSplitter()

    def initUI(self):
        title = "Convergence Viewer"
        if self.gene_name:
            title += f": {self.gene_name}"
        self.setWindowTitle(title)
        main_layout = QVBoxLayout(self)

        # Keep references to these checkboxes as None until needed
        self.onlyCcsCheck = None
        self.hideControlConvCheck = None

        # --- TOP CONTROLS ROW ---
        self.top_hbox = QHBoxLayout()  # Must be an attribute so we can reference it later

        font_bold = QFont()
        font_bold.setBold(True)

        # Label + Sort Combo
        sort_label = QLabel("Sort Sites By:")
        sort_label.setFont(font_bold)
        self.sort_combo = QComboBox()
        self.sort_combo.addItem("High Score → Low")   # idx 0
        self.sort_combo.addItem("Position Ascending") # idx 1
        self.top_hbox.addWidget(sort_label)
        self.top_hbox.addWidget(self.sort_combo)

        # "Show All Species" checkbox is kept for layout index calculations but hidden from the UI.
        self.showAllCheck = QCheckBox("Show All Species")
        self.showAllCheck.setChecked(True)  # always on
        self.showAllCheck.setVisible(False)  # hide from users
        # No stateChanged connection so users cannot toggle
        self.top_hbox.addWidget(self.showAllCheck)

        # --- Filters dropdown ---
        filters_label = QLabel("Filters:")
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("No Filter")                    # idx 0
        self.filter_combo.addItem("Show Only Selected Sites")     # idx 1
        self.filter_combo.addItem("Show Only CCS Sites")          # idx 2

        # Disable options based on data availability
        self.filter_combo.model().item(1).setEnabled(self.has_selected_sites)
        self.filter_combo.model().item(2).setEnabled(self.has_all_three)

        # Default selection
        default_index = 1 if self.has_selected_sites else 0
        self.filter_combo.setCurrentIndex(default_index)
        # Initialise flags BEFORE connecting signal
        self.only_selected = (default_index == 1)
        self.only_ccs      = (default_index == 2)

        self.filter_combo.currentIndexChanged.connect(self.onFilterChanged)

        self.top_hbox.addWidget(filters_label)
        self.top_hbox.addWidget(self.filter_combo)

        # Placeholders for backward compatibility (old radio buttons removed)
        self.onlySelectedCheck = None
        self.onlyCcsCheck = None

        # Hide Control Convergence checkbox (independent)
        self.hideControlConvCheck = QCheckBox("Hide Control Convergence")
        self.hideControlConvCheck.stateChanged.connect(self.onHideControlConvChanged)
        self.top_hbox.addWidget(self.hideControlConvCheck)
        # Set initial visibility based on current filter
        self.updateHideControlConvVisibility()

        # Helper to refresh enabled/disabled state based on data availability
        self.updateFilterComboStates()

        # Help button on far right
        self.help_button = QPushButton("Help")
        self.help_button.clicked.connect(self.showHelp)
        self.top_hbox.addStretch(1)
        self.top_hbox.addWidget(self.help_button)

        # ----------------------------------------------------------------------------------
        # END of UI setup -----------------------------------------------------------------
        # ----------------------------------------------------------------------------------

        main_layout.addLayout(self.top_hbox)

        # --- TABLES (TOP SPLITTER) ---
        self.top_left_table = QTableWidget()
        self.top_left_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.top_left_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.top_left_table.verticalHeader().setVisible(False)
        self.top_left_table.horizontalHeader().setVisible(False)
        self.top_left_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.top_left_table.customContextMenuRequested.connect(
            lambda pos: self.onSpeciesPaneContextMenu(pos, top=True)
        )
        self.top_left_table.horizontalHeader().setStretchLastSection(True)
        self.top_left_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        self.top_right_table = QTableWidget()
        self.top_right_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.top_right_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.top_right_table.verticalHeader().setVisible(False)
        self.top_right_table.horizontalHeader().setVisible(False)
        self.top_right_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        # Sync vertical scrolling
        self.top_left_table.verticalScrollBar().valueChanged.connect(
            self.top_right_table.verticalScrollBar().setValue
        )
        self.top_right_table.verticalScrollBar().valueChanged.connect(
            self.top_left_table.verticalScrollBar().setValue
        )
        # no horizontal scroll on the left
        self.top_left_table.horizontalScrollBar().setDisabled(True)

        self.top_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.top_splitter.addWidget(self.top_left_table)
        self.top_splitter.addWidget(self.top_right_table)

        # --- TABLES (BOTTOM SPLITTER) ---
        self.bottom_left_table = QTableWidget()
        self.bottom_left_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.bottom_left_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.bottom_left_table.verticalHeader().setVisible(False)
        self.bottom_left_table.horizontalHeader().setVisible(False)
        self.bottom_left_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.bottom_left_table.customContextMenuRequested.connect(
            lambda pos: self.onSpeciesPaneContextMenu(pos, top=False)
        )
        self.bottom_left_table.horizontalHeader().setStretchLastSection(True)
        self.bottom_left_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        self.bottom_right_table = QTableWidget()
        self.bottom_right_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.bottom_right_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.bottom_right_table.verticalHeader().setVisible(False)
        self.bottom_right_table.horizontalHeader().setVisible(False)
        self.bottom_right_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)

        # Sync vertical scrolling
        self.bottom_left_table.verticalScrollBar().valueChanged.connect(
            self.bottom_right_table.verticalScrollBar().setValue
        )
        self.bottom_right_table.verticalScrollBar().valueChanged.connect(
            self.bottom_left_table.verticalScrollBar().setValue
        )
        # Sync horizontal scrolling top->bottom
        self.top_right_table.horizontalScrollBar().valueChanged.connect(
            self.bottom_right_table.horizontalScrollBar().setValue
        )
        self.bottom_right_table.horizontalScrollBar().valueChanged.connect(
            self.top_right_table.horizontalScrollBar().setValue
        )

        self.bottom_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.bottom_splitter.addWidget(self.bottom_left_table)
        self.bottom_splitter.addWidget(self.bottom_right_table)

        self.vertical_splitter = QSplitter(Qt.Orientation.Vertical)
        # Let the bottom pane take extra vertical space on resize
        self.vertical_splitter.setStretchFactor(0, 0)  # top minimal but fixed proportion
        self.vertical_splitter.setStretchFactor(1, 1)  # bottom grows
        self.vertical_splitter.addWidget(self.top_splitter)
        self.vertical_splitter.addWidget(self.bottom_splitter)

        # if user unchecks 'Show All Species', we hide bottom
        self.bottom_splitter.setVisible(True)

        self.top_splitter.splitterMoved.connect(self.syncHorizontalSplitter)
        self.bottom_splitter.splitterMoved.connect(self.syncHorizontalSplitter)

        self.top_splitter.setSizes([240, 960])
        # Ensure bottom splitter starts at the same position as the top splitter
        self.bottom_splitter.setSizes(self.top_splitter.sizes())
        main_layout.addWidget(self.vertical_splitter, stretch=1)

        # --- THRESHOLD SLIDER + HISTOGRAM ---
        bottom_hbox = QHBoxLayout()
        slider_box = QVBoxLayout()

        slider_label = QLabel("Convergence Score Threshold:")
        slider_font = QFont()
        slider_font.setBold(True)
        slider_label.setFont(slider_font)

        self.threshold_slider = QSlider(Qt.Orientation.Horizontal)
        self.threshold_slider.setMinimum(self.min_score)
        self.threshold_slider.setMaximum(self.max_score)
        self.threshold_slider.setValue(self.default_threshold)
        self.threshold_slider.setFixedWidth(300)
        self.threshold_slider.valueChanged.connect(self.onThresholdChanged)

        self.current_threshold_label = QLabel(str(self.default_threshold))

        slider_box.addWidget(slider_label)
        slider_box.addWidget(self.threshold_slider)
        slider_box.addWidget(self.current_threshold_label)
        slider_box.addStretch(1)

        # 1) Create labels for protein length and sites shown
        self.protein_length_label = QLabel(f"Protein Length: {self.seq_length}")
        self.num_sites_shown_label = QLabel("Number of Sites Shown: 0")

        # 2) Add them to the same vertical layout with the threshold slider
        slider_box.addWidget(self.protein_length_label)
        slider_box.addWidget(self.num_sites_shown_label)

        self.hist_canvas = HistogramCanvas(self, width=4, height=2, dpi=100)
        self.hist_canvas.plot_scores(self.scores, self.default_threshold)

        bottom_hbox.addLayout(slider_box)
        bottom_hbox.addWidget(self.hist_canvas, stretch=1)
        main_layout.addLayout(bottom_hbox, stretch=0)

        self.setLayout(main_layout)
        # Provide a taller default window to better fit all tables
        self.resize(1200, 1000)

        # Default sort -> Position Ascending
        self.sort_combo.blockSignals(True)
        self.sort_combo.setCurrentIndex(1)
        self.default_sort_mode = "position"
        self.sort_combo.blockSignals(False)
        self.sort_combo.currentIndexChanged.connect(self.onSortModeChanged)


    def showHelp(self):
        """Help dialog content for the ESL-PSC Site Viewer."""
        msg = (
            "<h3>ESL-PSC Site Viewer Help</h3>"
            "<b>Purpose</b>: Explore per-site convergence statistics generated by an ESL-PSC run.<br><br>"
            "<b>Species Categories</b>: Convergent, Control, Outgroup and Other.<br>"
            "Right-click a species name to move it between Categories. Whenever Categories membership changes, the scores, CCS flags and histogram refresh automatically.<br><br>"
            "<b>Filters</b>:<br>"
            "• <i>No Filter</i> – show all sites.<br>"
            "• <i>Show Only Selected Sites</i> – show sites listed in the PSS (Predicted Selected Sites) file.<br>"
            "• <i>Show Only CCS Sites</i> – available only when all three focal Categories exist; shows sites meeting the CCS definition below.<br>"
            "• <i>Hide Control Convergence</i> – visible when all three focal Categories exist; hides sites where the Convergent residue matches the Outgroup residue and ≥2 Control species share a different residue.<br><br>"
            "<b>Convergence Score</b>: For each alignment column we sum the absolute differences in amino-acid usage between the Convergent and Control groups, then subtract 1 for every gap ('-') in those two groups (minimum 0). Gaps in the Outgroup are ignored for this penalty.<br><br>"
            "<b>CCS (Convergence at Conservative Sites)</b>: All Control and Outgroup species share residue R and at least two Convergent species share a residue ≠ R. Such sites are flagged with an asterisk. See: Xu et al. 2017, Mol. Biol. Evol., https://doi.org/10.1093/molbev/msw277<br>"
            "<b>Control-Convergence</b>: All Convergent match Outgroup residue R and at least two Control species share a residue ≠ R. Use the checkbox to hide these sites.<br><br>"
            "<b>Other Controls & Tips</b>:<br>"
            "• Sort Sites By – position (ascending) or score (descending).<br>"
            "• Threshold slider – show only sites with score ≥ chosen threshold (default 2; the red dashed line in the histogram marks the threshold).<br>"
            "• Show All Species – toggles visibility of the 'Other' species table.<br><br>"
            "See the ESL-PSC documentation for further details."
        )
        QMessageBox.information(self, "Help: ESL-PSC Site Viewer", msg)

    def syncHorizontalSplitter(self, pos, index):
        if self._syncing_horizontal_splitters:
            return
        self._syncing_horizontal_splitters = True
        sender_splitter = self.sender()
        if sender_splitter is self.top_splitter:
            self.bottom_splitter.setSizes(self.top_splitter.sizes())
        else:
            self.top_splitter.setSizes(self.bottom_splitter.sizes())
        self._syncing_horizontal_splitters = False

    def determine_default_threshold(self, num_sites_target=10):
        if not self.all_sites_info:
            return 0
        smin = min(self.scores)
        smax = max(self.scores)
        if smin == smax:
            return smin
        sorted_info = sorted(self.all_sites_info, key=lambda x: x['converge_degree'], reverse=True)
        if len(sorted_info) <= num_sites_target:
            return sorted_info[-1]['converge_degree']
        cscore_at_target = sorted_info[num_sites_target - 1]['converge_degree']
        i = num_sites_target
        while i < len(sorted_info):
            if sorted_info[i]['converge_degree'] == cscore_at_target:
                i += 1
            else:
                break
        return cscore_at_target

    def onSortModeChanged(self, idx):
        self.default_sort_mode = "score" if idx == 0 else "position"
        self.rebuildTables()

    def onShowAllChanged(self, state):
        """
        Toggles visibility of the 'other species' bottom splitter.
        """
        # ``state`` arrives as an ``int`` while ``Qt.CheckState.Checked`` is an
        # enum in PyQt6.  Cast explicitly so the comparison works on PyQt5/6.
        self.show_all_species = (Qt.CheckState(state) == Qt.CheckState.Checked)
        if self.show_all_species:
            self.bottom_splitter.setSizes(self.top_splitter.sizes())
        self.bottom_splitter.setVisible(True)
        self.rebuildTables()

    def onOnlyCcsChanged(self, checked):
        """Radio toggled for CCS filter."""
        self.only_ccs = bool(checked)
        # ensure mutual exclusivity handled by button group
        self.rebuildTables()

    def onHideControlConvChanged(self, state):
        # ``state`` is an ``int``; use explicit enum conversion for PyQt6
        self.hide_control_convergence = (
            Qt.CheckState(state) == Qt.CheckState.Checked
        )
        self.rebuildTables()

    def onOnlySelectedChanged(self, checked):
        """Radio toggled for selected sites filter."""
        self.only_selected = bool(checked)
        self.rebuildTables()

    # --------------------------- Filter helpers ---------------------------
    def onFilterChanged(self, idx):
        """Handle dropdown selection change."""
        self.only_selected = (idx == 1)
        self.only_ccs = (idx == 2)
        self.updateHideControlConvVisibility()
        self.rebuildTables()

    def updateHideControlConvVisibility(self):
        # Visible whenever all three species groups are present
        self.hideControlConvCheck.setVisible(self.has_all_three)
        # Set checked by default when visible
        if self.has_all_three:
            self.hideControlConvCheck.setChecked(True)
            self.hide_control_convergence = True

    def updateFilterComboStates(self):
        # Enable or disable dropdown items based on data availability
        self.filter_combo.model().item(1).setEnabled(self.has_selected_sites)
        self.filter_combo.model().item(2).setEnabled(self.has_all_three)
        # If current selection is disabled, reset to first valid option
        if not self.filter_combo.model().item(self.filter_combo.currentIndex()).isEnabled():
            new_idx = 1 if self.has_selected_sites else 0
            self.filter_combo.setCurrentIndex(new_idx)

    # --------------------------------------------------------------------

    def onThresholdChanged(self, val):
        self.current_threshold = val
        self.current_threshold_label.setText(str(val))
        self.hist_canvas.plot_scores(self.scores, val)
        self.rebuildTables()

    def onSpeciesPaneContextMenu(self, pos, top=True):
        """
        Right-click menu for species in the left table. 
        Dynamically adds or removes the 'Show Only CCS Sites' and 
        'Hide Control Convergence' checkboxes as groups 
        appear or disappear.
        """
        table = self.top_left_table if top else self.bottom_left_table
        row = table.rowAt(pos.y())
        if row < 0:
            return

        item = table.item(row, 0)
        if not item:
            return

        display_txt = item.text().strip()
        # Remove phenotype annotation such as ' (C4)' if present
        sp_name = display_txt.split(" (")[0].strip()
        # skip if label or blank
        if sp_name in ["Convergent Species", "Control Species", "Outgroup Species", ""]:
            return

        # figure out the current group
        if sp_name in self.convergent_species:
            group = "convergent"
        elif sp_name in self.control_species:
            group = "control"
        elif sp_name in self.outgroup_species:
            group = "outgroup"
        else:
            group = "other"

        menu = QMenu(self)
        if group == "convergent":
            toCtrl = menu.addAction("Move to Control")
            toOut = menu.addAction("Move to Outgroup")
            toOther = menu.addAction("Move to Other")
            chosen = menu.exec(table.mapToGlobal(pos))
            if chosen == toCtrl:
                self.convergent_species.remove(sp_name)
                self.control_species.append(sp_name)
                self.control_species.sort()
            elif chosen == toOut:
                self.convergent_species.remove(sp_name)
                self.outgroup_species.append(sp_name)
                self.outgroup_species.sort()
            elif chosen == toOther:
                self.convergent_species.remove(sp_name)

        elif group == "control":
            toConv = menu.addAction("Move to Convergent")
            toOut = menu.addAction("Move to Outgroup")
            toOther = menu.addAction("Move to Other")
            chosen = menu.exec(table.mapToGlobal(pos))
            if chosen == toConv:
                self.control_species.remove(sp_name)
                self.convergent_species.append(sp_name)
                self.convergent_species.sort()
            elif chosen == toOut:
                self.control_species.remove(sp_name)
                self.outgroup_species.append(sp_name)
                self.outgroup_species.sort()
            elif chosen == toOther:
                self.control_species.remove(sp_name)

        elif group == "outgroup":
            toConv = menu.addAction("Move to Convergent")
            toCtrl = menu.addAction("Move to Control")
            toOther = menu.addAction("Move to Other")
            chosen = menu.exec(table.mapToGlobal(pos))
            if chosen == toConv:
                self.outgroup_species.remove(sp_name)
                self.convergent_species.append(sp_name)
                self.convergent_species.sort()
            elif chosen == toCtrl:
                self.outgroup_species.remove(sp_name)
                self.control_species.append(sp_name)
                self.control_species.sort()
            elif chosen == toOther:
                self.outgroup_species.remove(sp_name)

        else:  # group == "other"
            toConv = menu.addAction("Move to Convergent")
            toCtrl = menu.addAction("Move to Control")
            toOut = menu.addAction("Move to Outgroup")
            chosen = menu.exec(table.mapToGlobal(pos))
            if chosen == toConv:
                self.convergent_species.append(sp_name)
                self.convergent_species.sort()
            elif chosen == toCtrl:
                self.control_species.append(sp_name)
                self.control_species.sort()
            elif chosen == toOut:
                self.outgroup_species.append(sp_name)
                self.outgroup_species.sort()

        # check if we have all 3 now
        old_has_all_three = self.has_all_three
        new_has_all_three = (
            bool(self.convergent_species)
            and bool(self.control_species)
            and bool(self.outgroup_species)
        )
        self.has_all_three = new_has_all_three

        if new_has_all_three and not old_has_all_three:
            # CCS filter becomes available
            self.updateFilterComboStates()
            self.updateHideControlConvVisibility()
            # All UI widgets already exist; just update dropdown states

        elif not new_has_all_three and old_has_all_three:
            # Lost a group; CCS filter no longer available
            self.updateFilterComboStates()
            self.updateHideControlConvVisibility()
            # No need to remove any widgets; just reset states as needed
            self.only_ccs = False
            self.hide_control_convergence = False

        # Update remembered outgroup list after any change
        SiteViewer.REMEMBERED_OUTGROUP = self.outgroup_species.copy()

        # now recalc scores + rebuild
        self.recalc_scores()
        self.rebuildTables()

    # ------------------------------------------------------------------
    # Ensure the remembered outgroup list is persisted when the window closes
    # ------------------------------------------------------------------
    def closeEvent(self, event):  # noqa: N802 (Qt override)
        SiteViewer.REMEMBERED_OUTGROUP = self.outgroup_species.copy()
        super().closeEvent(event)

    def recalc_scores(self):
        """
        Recalculate the convergence scores for each site.
        Singletons (non-gap) in Convergent+Control are replaced by '?'.
        Each gap ('-') in Convergent or Control reduces the score by 1.
        """
        # Initialize threshold_slider if it doesn't exist yet
        if not hasattr(self, 'threshold_slider'):
            return  # UI not fully initialized yet

        conv_indices = [self.species_ids.index(sp) for sp in self.convergent_species if sp in self.species_ids]
        ctrl_indices = [self.species_ids.index(sp) for sp in self.control_species if sp in self.species_ids]
        out_indices  = [self.species_ids.index(sp) for sp in self.outgroup_species if sp in self.species_ids]

        updated_sites = []

        for old_info in self.all_sites_info:
            pos = old_info['position']

            residues = [seq[pos] for seq in self.sequences]

            # Identify singletons among Convergent+Control ignoring gaps
            cc_indices = conv_indices + ctrl_indices
            cc_residues = [residues[i] for i in cc_indices]
            cc_counted = Counter(cc_residues)

            no_singletons = residues[:]
            for i in cc_indices:
                if cc_counted[residues[i]] == 1 and residues[i] != '-':
                    no_singletons[i] = '?'

            conv_aa = [no_singletons[i] for i in conv_indices]
            ctrl_aa = [no_singletons[i] for i in ctrl_indices]
            out_aa  = [no_singletons[i] for i in out_indices]

            # Exclude both '?' and '-' from raw difference calculation
            conv_clean = [r for r in conv_aa if r not in ('?', '-')]
            ctrl_clean = [r for r in ctrl_aa if r not in ('?', '-')]
            ccounts = Counter(conv_clean)
            ccounts.subtract(ctrl_clean)
            raw_score = sum(abs(v) for v in ccounts.values())

            # Gap penalty: count each '-' in conv+ctrl and subtract 1 per gap
            gap_count = sum(r == '-' for r in conv_aa) + sum(r == '-' for r in ctrl_aa)
            final_score = raw_score - gap_count
            if final_score < 0:
                final_score = 0

            # CCS / Control-Convergence detection if we have all 3 groups
            is_ccs = False
            conv_same_as_ctrl = False
            if self.has_all_three:
                # Cleaned residue lists (exclude gaps and '?')
                clean_conv = [x for x in conv_aa if x not in ('?', '-')]
                clean_ctrl = [x for x in ctrl_aa if x not in ('?', '-')]
                clean_out  = [x for x in out_aa  if x not in ('?', '-')]

                # --------------------------------------------------
                # CCS site (original definition):
                #   Control & Outgroup share residue R
                #   ≥2 Convergent share residue != R
                # --------------------------------------------------
                if (
                    clean_ctrl and clean_out and
                    len(set(clean_ctrl)) == 1 and len(set(clean_out)) == 1 and
                    list(set(clean_ctrl))[0] == list(set(clean_out))[0]
                ):
                    ctrl_res = clean_ctrl[0]
                    conv_counter = Counter(clean_conv)
                    for res, cnt in conv_counter.items():
                        if res != ctrl_res and cnt >= 2:
                            is_ccs = True
                            break

                # --------------------------------------------------
                # Control-convergence (inverse of CCS):
                #   Convergent group matches Outgroup residue R
                #   ≥2 Control species share residue ≠ R
                # --------------------------------------------------
                if clean_out and len(set(clean_out)) == 1:
                    out_res = clean_out[0]
                    # All convergent match out_res
                    if clean_conv and all(r == out_res for r in clean_conv):
                        ctrl_counter = Counter(clean_ctrl)
                        for res, cnt in ctrl_counter.items():
                            if res != out_res and cnt >= 2:
                                conv_same_as_ctrl = True
                                break


            new_info = {
                'position': pos,
                'converge_degree': final_score,
                'is_ccs': is_ccs,
                'conv_same_as_ctrl': conv_same_as_ctrl
            }
            updated_sites.append(new_info)

        self.all_sites_info = updated_sites
        self.scores = [s['converge_degree'] for s in updated_sites]

        # Update slider bounds to reflect new score range
        self.min_score = min(self.scores) if self.scores else 0
        self.max_score = max(self.scores) if self.scores else 0
        self.threshold_slider.setMinimum(self.min_score)
        self.threshold_slider.setMaximum(self.max_score)

        # Auto-align threshold once at start
        if not self._initial_threshold_aligned:
            if not (self.min_score <= self.current_threshold <= self.max_score):
                self.current_threshold = max(self.min_score, min(self.max_score, self.default_threshold))
            self.threshold_slider.blockSignals(True)
            self.threshold_slider.setValue(self.current_threshold)
            self.threshold_slider.blockSignals(False)
            self.current_threshold_label.setText(str(self.current_threshold))
            self._initial_threshold_aligned = True
        else:
            # Keep threshold within bounds if scores changed drastically
            if self.current_threshold < self.min_score:
                self.current_threshold = self.min_score
                self.threshold_slider.setValue(self.current_threshold)
            elif self.current_threshold > self.max_score:
                self.current_threshold = self.max_score
                self.threshold_slider.setValue(self.current_threshold)

        # Redraw histogram using the (possibly updated) threshold
        self.hist_canvas.plot_scores(self.scores, self.current_threshold)
        

    def rebuildTables(self):
        # Skip if UI isn't fully initialized yet
        if not hasattr(self, 'threshold_slider'):
            return
            
        # recalc with updated gap penalty and ccs logic
        self.recalc_scores()

        # apply threshold, ccs/hide filters
        filtered_sites = []
        for s in self.all_sites_info:
            if s['converge_degree'] < self.current_threshold:
                continue
            if self.has_all_three and self.only_ccs and not s['is_ccs']:
                continue
            if self.has_all_three and self.hide_control_convergence and s['conv_same_as_ctrl']:
                continue
            if self.only_selected and s['position'] not in self.pss_scores:
                continue
            filtered_sites.append(s)

        if self.default_sort_mode == "score":
            filtered_sites.sort(key=lambda x: (-x['converge_degree'], x['position']))
        else:
            filtered_sites.sort(key=lambda x: x['position'])

        # Update the "Number of Sites Shown" label
        self.num_sites_shown_label.setText(f"Number of Sites Shown: {len(filtered_sites)}")


        self._buildTopTables(filtered_sites)
        self._buildBottomTables(filtered_sites)
        self._unifyCols()
        # Re-adjust vertical splitter after tables rebuild
        self._adjustVerticalSplitter()

    def _buildTopTables(self, displayed_sites):
        self.top_left_table.clearContents()
        self.top_right_table.clearContents()

        conv_size = len(self.convergent_species)
        ctrl_size = len(self.control_species)
        out_size = len(self.outgroup_species)

        row_idx = 0
        # row0 => "Position", row1 => "Score", optional row2 => "PSS", then blank
        base_rows = 3 if not self.pss_scores else 4
        row_idx += base_rows

        conv_label_row = None
        if conv_size > 0:
            conv_label_row = row_idx
            row_idx += 1 + conv_size

        if conv_label_row:
            row_idx += 1

        ctrl_label_row = None
        if ctrl_size > 0:
            ctrl_label_row = row_idx
            row_idx += 1 + ctrl_size

        if ctrl_label_row:
            row_idx += 1

        out_label_row = None
        if out_size > 0:
            out_label_row = row_idx
            row_idx += 1 + out_size

        total_rows = row_idx
        if total_rows < 3:
            total_rows = 3

        self.top_left_table.setRowCount(total_rows)
        self.top_left_table.setColumnCount(1)

        bold_f = QFont()
        bold_f.setBold(True)

        def left_item(txt, b=False):
            i = QTableWidgetItem(txt)
            i.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            if b:
                i.setFont(bold_f)
            return i

        # row0 => Position
        self.top_left_table.setItem(0, 0, left_item("Position", True))
        # row1 => Score
        self.top_left_table.setItem(1, 0, left_item("Convergence Score", True))
        next_row = 2
        if self.pss_scores:
            self.top_left_table.setItem(2, 0, left_item("PSS", True))
            next_row = 3
        # blank separator
        self.top_left_table.setItem(next_row, 0, left_item(""))

        # conv
        if conv_label_row is not None:
            self.top_left_table.setItem(conv_label_row, 0, left_item("Convergent Species", True))
            for i, sp in enumerate(self.convergent_species):
                rr = conv_label_row + 1 + i
                self.top_left_table.setItem(rr,0, self._create_species_item(sp))

        # ctrl
        if ctrl_label_row is not None:
            self.top_left_table.setItem(ctrl_label_row, 0, left_item("Control Species", True))
            for i, sp in enumerate(self.control_species):
                rr = ctrl_label_row + 1 + i
                self.top_left_table.setItem(rr,0, self._create_species_item(sp))

        # outgroup
        if out_label_row is not None:
            self.top_left_table.setItem(out_label_row, 0, left_item("Outgroup Species", True))
            for i, sp in enumerate(self.outgroup_species):
                rr = out_label_row + 1 + i
                self.top_left_table.setItem(rr,0, self._create_species_item(sp))

        for r in range(total_rows):
            self.top_left_table.setRowHeight(r, 24)

        self.top_right_table.setRowCount(total_rows)
        displayed_cols = len(displayed_sites)
        self.top_right_table.setColumnCount(displayed_cols)

        # row0 => position, row1 => score, optional PSS row, then blank
        for c, site in enumerate(displayed_sites):
            pos_str = str(site['position']+1)
            pi = QTableWidgetItem(pos_str)
            pi.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.top_right_table.setItem(0,c, pi)

            sc_val = site['converge_degree']
            sc_str = f"{sc_val}*" if site['is_ccs'] else str(sc_val)
            si = QTableWidgetItem(sc_str)
            si.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.top_right_table.setItem(1, c, si)

            if self.pss_scores:
                pss_val = self.pss_scores.get(site['position'])
                pss_str = f"{pss_val:.3f}" if pss_val is not None else ""
                pi2 = QTableWidgetItem(pss_str)
                pi2.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.top_right_table.setItem(2, c, pi2)

        # conv
        if conv_label_row is not None:
            for i, sp in enumerate(self.convergent_species):
                rr = conv_label_row + 1 + i
                seq = self.get_seq_for_species(sp)
                if seq:
                    for c, site in enumerate(displayed_sites):
                        aa = seq[site['position']]
                        it = self.make_aa_item(aa)
                        self.top_right_table.setItem(rr,c, it)

        # ctrl
        if ctrl_label_row is not None:
            for i, sp in enumerate(self.control_species):
                rr = ctrl_label_row + 1 + i
                seq = self.get_seq_for_species(sp)
                if seq:
                    for c, site in enumerate(displayed_sites):
                        aa = seq[site['position']]
                        it = self.make_aa_item(aa)
                        self.top_right_table.setItem(rr,c, it)

        # outgroup
        if out_label_row is not None:
            for i, sp in enumerate(self.outgroup_species):
                rr = out_label_row + 1 + i
                seq = self.get_seq_for_species(sp)
                if seq:
                    for c, site in enumerate(displayed_sites):
                        aa = seq[site['position']]
                        it = self.make_aa_item(aa)
                        self.top_right_table.setItem(rr,c, it)

        for r in range(total_rows):
            self.top_right_table.setRowHeight(r, 24)

    def _buildBottomTables(self, displayed_sites):
        # "other" species
        used_set = set(self.convergent_species + self.control_species + self.outgroup_species)
        other_species = [sp for sp in self.all_species if sp not in used_set]
        if not other_species:
            self.bottom_left_table.setRowCount(0)
            self.bottom_right_table.setRowCount(0)
            return

        self.bottom_left_table.clearContents()
        self.bottom_right_table.clearContents()

        bold_f = QFont()
        bold_f.setBold(True)

        row0_label = 0
        row_sp_start = 1
        row_sp_end = row_sp_start + len(other_species)-1
        total_rows = row_sp_end+1
        if total_rows < 1:
            total_rows = 1

        self.bottom_left_table.setRowCount(total_rows)
        self.bottom_left_table.setColumnCount(1)

        def left_item(txt, b=False):
            i = QTableWidgetItem(txt)
            i.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            if b:
                i.setFont(bold_f)
            return i

        def species_label_helper(sp):
            disp_txt = sp
            if sp in self.species_pheno_map:
                p_val = self.species_pheno_map[sp]
                p_name = self.pheno_name_map.get(p_val, str(p_val))
                disp_txt += f" ({p_name})"
            item = left_item(disp_txt)
            if sp in self.species_pheno_map and self.species_pheno_map[sp] == self.convergent_pheno_value:
                item.setForeground(QBrush(QColor("blue")))
            return item

        self.bottom_left_table.setItem(row0_label, 0, left_item("Other Species", True))

        for i, sp in enumerate(other_species):
            rr = row_sp_start + i
            self.bottom_left_table.setItem(rr, 0, self._create_species_item(sp))

        for r in range(total_rows):
            self.bottom_left_table.setRowHeight(r, 24)

        self.bottom_right_table.setRowCount(total_rows)
        cols = len(displayed_sites)
        self.bottom_right_table.setColumnCount(cols)

        for r in range(total_rows):
            for c in range(cols):
                self.bottom_right_table.setItem(r,c, QTableWidgetItem(""))

        for c in range(cols):
            self.bottom_right_table.setItem(row0_label,c, QTableWidgetItem(""))

        for i, sp in enumerate(other_species):
            rr = row_sp_start + i
            seq = self.get_seq_for_species(sp)
            if seq:
                for cc, site in enumerate(displayed_sites):
                    aa = seq[site['position']]
                    it = self.make_aa_item(aa)
                    self.bottom_right_table.setItem(rr,cc, it)

        for r in range(total_rows):
            self.bottom_right_table.setRowHeight(r,24)

    def _create_species_item(self, sp: str) -> QTableWidgetItem:
        """Return a table item for a species with phenotype annotation and coloring."""
        disp_txt = sp
        if sp in self.species_pheno_map:
            p_val = self.species_pheno_map[sp]
            p_name = self.pheno_name_map.get(p_val, str(p_val))
            disp_txt += f" ({p_name})"
        item = QTableWidgetItem(disp_txt)
        item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        if sp in self.species_pheno_map and self.species_pheno_map[sp] == self.convergent_pheno_value:
            item.setForeground(QBrush(QColor("blue")))
        return item

    def _adjustVerticalSplitter(self):
        """Set vertical splitter sizes so top pane shows conv+ctrl species without scroll if space allows."""
        if not self.top_left_table.rowCount():
            return
        rows = self.top_left_table.rowCount()
        row_h = self.top_left_table.rowHeight(0) or 24
        top_needed = rows * row_h + 40  # padding for headers etc.
        bottom_min = 80  # small minimum height for other species pane
        total_height = max(self.vertical_splitter.height(), top_needed + bottom_min)
        # Ensure top isn't bigger than total - bottom_min
        top_size = min(top_needed, total_height - bottom_min)
        bottom_size = total_height - top_size
        self.vertical_splitter.setSizes([top_size, bottom_size])


    def _unifyCols(self):
        self.top_right_table.resizeColumnsToContents()
        tcols = self.top_right_table.columnCount()
        twidths = [self.top_right_table.columnWidth(c) for c in range(tcols)]

        self.bottom_right_table.resizeColumnsToContents()
        bcols = self.bottom_right_table.columnCount()
        bwidths = [self.bottom_right_table.columnWidth(c) for c in range(bcols)]

        combined = twidths + bwidths
        if not combined:
            return
        max_w = max(combined)

        for c in range(tcols):
            self.top_right_table.setColumnWidth(c, max_w)
        for c in range(bcols):
            self.bottom_right_table.setColumnWidth(c, max_w)

    def make_aa_item(self, aa):
        it = QTableWidgetItem(aa)
        it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        color_hex = ZAPPO_STATIC_COLORS.get(aa.upper(), "#C8C8C8")
        color = QColor(color_hex)
        avg_rgb = (color.red()+color.green()+color.blue())/3
        text_c = QColor(0,0,0) if avg_rgb>=128 else QColor(255,255,255)
        it.setBackground(QBrush(color))
        it.setForeground(QBrush(text_c))
        return it

    def get_seq_for_species(self, sid):
        for r_sid, seq in self.records:
            if r_sid == sid:
                return seq
        return None
