# widgets/site_viewer.py
from __future__ import annotations

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush, QFont
from PyQt6.QtWidgets import (
    QTableWidget, QTableWidgetItem, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QSlider, QComboBox, QCheckBox, QLabel, QPushButton,
    QAbstractItemView, QMessageBox, QMenu
)

# point to your shared constants and canvas
from gui.ui.widgets.histogram_canvas import HistogramCanvas
from gui.constants import ZAPPO_STATIC_COLORS



class SiteViewer(QWidget):
    """
    Main alignment-inspection widget.
    """

    def __init__(
        self,
        records,                # List[tuple[id, seq]]
        convergent_species: List[str],
        control_species:    List[str],
        outgroup_species:   List[str],
        all_sites_info:     List[Dict[str, Any]] | None = None,
        show_all_by_default: bool = False,
        pss_scores: Dict[int, float] | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)

        # store the minimal bits up front
        self.records              = records
        self.convergent_species  = sorted(convergent_species)
        self.control_species     = sorted(control_species)
        self.outgroup_species    = sorted(outgroup_species)
        self.pss_scores          = pss_scores or {}

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

        self.default_threshold = self.determine_default_threshold()
        self.current_threshold = self.default_threshold

        self.default_sort_mode = "position"
        # Always show all species; the checkbox will be hidden but kept in the layout
        self.show_all_species = True

        # Additional controls if we have all three groups
        self.only_ccs = False
        self.hide_control_convergence = False

        # Are all 3 groups non-empty?
        self.has_all_three = bool(self.convergent_species) and bool(self.control_species) and bool(self.outgroup_species)

        self._syncing_horizontal_splitters = False

        self.initUI()
        self.rebuildTables()
        # Adjust vertical splitter sizes based on table content
        self._adjustVerticalSplitter()

    def initUI(self):
        self.setWindowTitle("Convergence Viewer")
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

        # If we start with all three groups, create CCS checkboxes right away
        if self.has_all_three:
            self.onlyCcsCheck = QCheckBox("Show Only CCS Sites")
            self.onlyCcsCheck.stateChanged.connect(self.onOnlyCcsChanged)
            self.top_hbox.addWidget(self.onlyCcsCheck)

            self.hideControlConvCheck = QCheckBox("Hide Control Convergence")
            self.hideControlConvCheck.stateChanged.connect(self.onHideControlConvChanged)
            self.top_hbox.addWidget(self.hideControlConvCheck)

        # Help button on far right
        self.help_button = QPushButton("Help")
        self.help_button.clicked.connect(self.showHelp)
        self.top_hbox.addStretch(1)
        self.top_hbox.addWidget(self.help_button)

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
        """Restored explanation of how the convergence score is calculated,
        including new gap penalty and other details.
        """
        msg = (
            "<b>Convergence Score Calculation</b><br>"
            "We examine alignment columns for differences between Convergent vs. Control usage, ignoring singletons. "
            "Each gap in the Convergent or Control group reduces the final score by 1 (minimum 0). "
            "Outgroup gaps are ignored for the gap penalty. The sum of absolute differences in amino-acid usage "
            "gives the raw score, and then we apply the gap penalty.<br><br>"
            "<b>CCS Sites</b>: If all Control share the same residue as all Outgroup, and at least 2 Convergent share "
            "a different residue, we mark that site with an asterisk in the score cell. CCS stands for \"Convergence at Conservative Sites\" "
            "(See: Xu, Shaohua, Ziwen He, Zixiao Guo, Zhang Zhang, Gerald J. Wyckoff, Anthony Greenberg, Chung-I Wu, and Suhua Shi. 2017."
            "“Genome-Wide Convergence during Evolution of Mangroves from Woody Plants.” Molecular Biology and Evolution 34 (4): 1008–15.)<br><br>"
            "<b>Other Controls</b>:<br>"
            "• Sort Sites By: Ascending position (default) or descending Convergence Score.<br>"
            "• Threshold Slider: Only show sites >= that score.<br>"
            "• Show All Species: Toggle the lower pane for 'other' species not in Convergent, Control, or Outgroup.<br>"
            "• Show Only CCS Sites: If all three species groups exist, this hides non-CCS sites.<br>"
            "• Hide Control Convergence: If all three species groups exist, hides sites where Convergent's dominant residue "
            "is the same as Control/Outgroup's residue.<br><br>"
            "• Right-click on the species names to move them between groups. "
            "When group membership changes, scores are recalculated."
        )
        QMessageBox.information(self, "Help: Convergence Viewer v1.1", msg)

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
        self.show_all_species = (state == Qt.CheckState.Checked)
        if self.show_all_species:
            self.bottom_splitter.setSizes(self.top_splitter.sizes())
        self.bottom_splitter.setVisible(True)
        self.rebuildTables()

    def onOnlyCcsChanged(self, state):
        self.only_ccs = (state == Qt.CheckState.Checked)
        self.rebuildTables()

    def onHideControlConvChanged(self, state):
        self.hide_control_convergence = (state == Qt.CheckState.Checked)
        self.rebuildTables()

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

        sp_name = item.text().strip()
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
            chosen = menu.exec_(table.mapToGlobal(pos))
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
            chosen = menu.exec_(table.mapToGlobal(pos))
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
            chosen = menu.exec_(table.mapToGlobal(pos))
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
            chosen = menu.exec_(table.mapToGlobal(pos))
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
            # Calculate the index of the "Show All Species" checkbox.
            # We want to insert new checkboxes right after this widget.
            show_all_index = self.top_hbox.indexOf(self.showAllCheck)
            insert_index = show_all_index + 1  # position immediately after "Show All Species"

            # Add "Show Only CCS Sites" checkbox if not already present
            if self.onlyCcsCheck is None:
                self.onlyCcsCheck = QCheckBox("Show Only CCS Sites")
                self.onlyCcsCheck.stateChanged.connect(self.onOnlyCcsChanged)
                # Insert the checkbox at our calculated position
                self.top_hbox.insertWidget(insert_index, self.onlyCcsCheck)
                insert_index += 1  # Increment index for next insertion

            # Add "Hide Control Convergence" checkbox if not already present
            if self.hideControlConvCheck is None:
                self.hideControlConvCheck = QCheckBox("Hide Control Convergence")
                self.hideControlConvCheck.stateChanged.connect(self.onHideControlConvChanged)
                # Insert next to the previous widget (still right after Show All Species)
                self.top_hbox.insertWidget(insert_index, self.hideControlConvCheck)

        elif not new_has_all_three and old_has_all_three:
            # We lost a group => remove checkboxes if they exist
            if self.onlyCcsCheck is not None:
                # Reset the internal 'only_ccs' state to False
                self.only_ccs = False

                self.top_hbox.removeWidget(self.onlyCcsCheck)
                self.onlyCcsCheck.deleteLater()
                self.onlyCcsCheck = None

            if self.hideControlConvCheck is not None:
                # Reset the internal 'hide_control_convergence' state to False
                self.hide_control_convergence = False

                self.top_hbox.removeWidget(self.hideControlConvCheck)
                self.hideControlConvCheck.deleteLater()
                self.hideControlConvCheck = None

        # now recalc scores + rebuild
        self.recalc_scores()
        self.rebuildTables()

    def recalc_scores(self):
        """
        Recalculate the convergence scores for each site.
        Singletons (non-gap) in Convergent+Control are replaced by '?'.
        Each gap ('-') in Convergent or Control reduces the score by 1.
        """
        from collections import Counter

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

            # CCS detection if we have all 3 groups
            is_ccs = False
            conv_same_as_ctrl = False
            if self.has_all_three:
                clean_ctrl = [x for x in ctrl_aa if x not in ('?', '-')]
                clean_out  = [x for x in out_aa if x not in ('?', '-')]

                if len(clean_ctrl) == len(ctrl_aa) and len(clean_out) == len(out_aa):
                    if len(set(clean_ctrl)) == 1 and len(set(clean_out)) == 1:
                        if list(set(clean_ctrl))[0] == list(set(clean_out))[0]:
                            clean_conv = [x for x in conv_aa if x not in ('?', '-')]
                            ctrl_res = clean_ctrl[0]
                            from collections import Counter
                            conv_counter = Counter(clean_conv)
                            for rrr, ccount in conv_counter.items():
                                if rrr != ctrl_res and ccount >= 2:
                                    is_ccs = True
                                    break

                # Hide-control-convergence detection:
                # Hide if all non-gap Convergent species match the shared Control/Outgroup residue and no singletons.
                if len(set(clean_ctrl)) == 1 and len(set(clean_out)) == 1:
                    shared_res = list(set(clean_ctrl))[0]
                    if shared_res == list(set(clean_out))[0]:
                        # Check that all non-gap species in Convergent match shared_res and there are no '?'
                        if all((x == shared_res or x == '-') for x in conv_aa) and '?' not in conv_aa:
                            conv_same_as_ctrl = True


            new_info = {
                'position': pos,
                'converge_degree': final_score,
                'is_ccs': is_ccs,
                'conv_same_as_ctrl': conv_same_as_ctrl
            }
            updated_sites.append(new_info)

        self.all_sites_info = updated_sites
        self.scores = [s['converge_degree'] for s in updated_sites]
        self.hist_canvas.plot_scores(self.scores, self.current_threshold)

        # Update slider range based on new scores
        self.min_score = min(self.scores) if self.scores else 0
        self.max_score = max(self.scores) if self.scores else 0
        self.threshold_slider.setMinimum(self.min_score)
        self.threshold_slider.setMaximum(self.max_score)

        # Ensure current threshold is within new bounds
        if self.current_threshold < self.min_score:
            self.current_threshold = self.min_score
            self.threshold_slider.setValue(self.current_threshold)
        elif self.current_threshold > self.max_score:
            self.current_threshold = self.max_score
            self.threshold_slider.setValue(self.current_threshold)
        

    def rebuildTables(self):
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

        blank_after_conv = row_idx
        if conv_label_row:
            row_idx += 1

        ctrl_label_row = None
        if ctrl_size > 0:
            ctrl_label_row = row_idx
            row_idx += 1 + ctrl_size

        blank_after_ctrl = row_idx
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
                self.top_left_table.setItem(rr,0, left_item(sp))

        # ctrl
        if ctrl_label_row is not None:
            self.top_left_table.setItem(ctrl_label_row, 0, left_item("Control Species", True))
            for i, sp in enumerate(self.control_species):
                rr = ctrl_label_row + 1 + i
                self.top_left_table.setItem(rr,0, left_item(sp))

        # outgroup
        if out_label_row is not None:
            self.top_left_table.setItem(out_label_row, 0, left_item("Outgroup Species", True))
            for i, sp in enumerate(self.outgroup_species):
                rr = out_label_row + 1 + i
                self.top_left_table.setItem(rr,0, left_item(sp))

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
                pss_str = f"{pss_val:.5f}" if pss_val is not None else ""
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

        self.bottom_left_table.setItem(row0_label,0, left_item("Other Species", True))

        for i, sp in enumerate(other_species):
            rr = row_sp_start + i
            self.bottom_left_table.setItem(rr,0, left_item(sp))

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
