from __future__ import annotations

from typing import Callable, Sequence

from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QDialogButtonBox,
    QLabel,
    QDoubleSpinBox,
    QPushButton,
    QSizePolicy,
    QComboBox,
    QCheckBox,
    QSpinBox,
    QRadioButton,
    QButtonGroup,
    QFileDialog,
    QLineEdit,
)

from gui.ui.widgets.histogram_canvas import HistogramCanvas


class AutoSelectOptionsDialog(QDialog):
    """Dialog to choose auto-select behavior with vertically stacked buttons.

    Presents options in this order: Longest Sequence → Shortest distance → Max trait contrast → Composite best → Random → Default → Cancel.
    Stores the chosen option as `self.choice` in
    {"default","longest","shortest","contrast","composite","random","pct_contrast",None}.
    """

    def __init__(
        self,
        allow_longest: bool,
        allow_contrast: bool,
        allow_composite: bool = True,
        allow_pct_contrast: bool = False,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Auto Select Options")
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Choose how to resolve equally valid siblings"))

        # --- Alternates controls ---
        # Number of alternates per convergent/control species (0..20)
        alts_row = QHBoxLayout()
        alts_row.addWidget(QLabel("Number of alternates to add per convergent and control species:"))
        self._alts_spin = QSpinBox()
        self._alts_spin.setRange(0, 20)
        self._alts_spin.setValue(0)
        self._alts_spin.setToolTip(
            "Select alternate siblings species in each pair if possible to allow multiple species combinations"
        )
        alts_row.addWidget(self._alts_spin)
        alts_row.addStretch()
        layout.addLayout(alts_row)

        # Maximum combinations (1..32,999). Disabled unless alternates > 0
        max_row = QHBoxLayout()
        max_lbl = QLabel("Maximum species combinations:")
        max_row.addWidget(max_lbl)
        self._max_spin = QSpinBox()
        self._max_spin.setRange(1, 32999)  # strictly less than 33,000
        self._max_spin.setValue(1)
        self._max_spin.setEnabled(False)
        self._max_spin.setToolTip("Upper limit on the product of choices across pairs. Disabled unless alternates > 0.")
        max_row.addWidget(self._max_spin)
        max_row.addStretch()
        layout.addLayout(max_row)

        # Enable/disable logic for max combinations
        def _on_alts_changed(val: int):
            try:
                if val <= 0:
                    self._max_spin.setValue(1)
                    self._max_spin.setEnabled(False)
                else:
                    # Enable and set a sensible default if previously disabled or at 1
                    was_enabled = self._max_spin.isEnabled()
                    self._max_spin.setEnabled(True)
                    if (not was_enabled) or self._max_spin.value() <= 1:
                        self._max_spin.setValue(32)
            except Exception:
                pass
        self._alts_spin.valueChanged.connect(_on_alts_changed)

        def mk_btn(text: str, choice: str, enabled: bool = True, tip: str | None = None) -> QPushButton:
            btn = QPushButton(text)
            btn.setEnabled(enabled)
            if tip:
                btn.setToolTip(tip)
            btn.clicked.connect(lambda: self._pick(choice))
            layout.addWidget(btn)
            return btn

        self.choice: str | None = None

        # Create buttons stacked vertically
        self.longest_btn = mk_btn(
            "Longest Sequence",
            "longest",
            enabled=True,
            tip=(
                "Choose tips with the greatest total sequence length under each ancestor. "
                "If no alignments directory is set, you will be prompted to choose one."
            ),
        )
        self.shortest_btn = mk_btn(
            "Shortest distance",
            "shortest",
            tip=(
                "Resolve equally valid siblings using a shortest-distance heuristic."
            ),
        )
        self.contrast_btn = mk_btn(
            "Max trait contrast",
            "contrast",
            enabled=allow_contrast,
            tip=(
                "Choose tips with maximal continuous-phenotype difference."
                if allow_contrast
                else "Available only when continuous phenotypes are loaded"
            ),
        )
        self.pct_contrast_btn = mk_btn(
            "Min % contrast (positive traits)",
            "pct_contrast",
            enabled=allow_pct_contrast,
            tip=(
                "Continuous positive traits only. Select pairs by closest distance first, "
                "requiring a minimum percent difference between upper and lower species."
                if allow_pct_contrast
                else "Available only when continuous phenotypes are loaded"
            ),
        )
        self.composite_btn = mk_btn(
            "Composite best",
            "composite",
            enabled=allow_composite,
            tip=(
                "Best pair of siblings based on a composite score of distance, and sequence length, "
                "and trait contrast if available (all rescaled)."
            ),
        )
        self.random_btn = mk_btn(
            "Random",
            "random",
            tip=("Choose randomly among siblings."),
        )
        self.default_btn = mk_btn(
            "Default",
            "default",
            tip=(
                "Use initial adjacent tips with opposite phenotypes."
            ),
        )

        # Cancel at the bottom
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        layout.addWidget(cancel)

    @property
    def num_alternates(self) -> int:
        try:
            return int(self._alts_spin.value())
        except Exception:
            return 0

    @property
    def max_combinations(self) -> int:
        try:
            return int(self._max_spin.value())
        except Exception:
            return 1

    def _pick(self, choice: str) -> None:
        self.choice = choice
        self.accept()


class PhenoThresholdDialog(QDialog):
    """Dialog for choosing thresholds for continuous phenotypes."""

    def __init__(self, values: Sequence[float], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phenotype Thresholds")
        self.setSizeGripEnabled(True)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        form = QHBoxLayout()
        form.addWidget(QLabel("Lower threshold:"))
        self.lower_spin = QDoubleSpinBox()
        form.addWidget(self.lower_spin)
        form.addWidget(QLabel("Upper threshold:"))
        self.upper_spin = QDoubleSpinBox()
        form.addWidget(self.upper_spin)
        layout.addLayout(form)

        # Quantile tail auto-set controls (manual numeric entry)
        qrow = QHBoxLayout()
        qlbl = QLabel("Auto-set by tails (%):")
        qlbl.setToolTip(
            "Enter a tail percentage to automatically set lower/upper thresholds to the corresponding\n"
            "quantiles (e.g., 10% sets thresholds at the 10th and 90th percentiles).\n"
            "Set to 0% to keep full manual control. At 50%, both thresholds are set to the median."
        )
        self.quantile_spin = QDoubleSpinBox()
        self.quantile_spin.setSuffix("%")
        self.quantile_spin.setRange(0.0, 50.0)  # cannot exceed 50%
        self.quantile_spin.setDecimals(1)
        self.quantile_spin.setSingleStep(1.0)
        self.quantile_spin.setToolTip(
            "0% keeps manual thresholds; >0% auto-applies symmetric quantile tails; 50% sets both to median."
        )
        self.quantile_spin.valueChanged.connect(self._on_quantile_changed)
        qrow.addWidget(qlbl)
        qrow.addWidget(self.quantile_spin)
        qrow.addStretch()
        layout.addLayout(qrow)

        self.values = list(values)
        if self.values:
            vmin, vmax = min(self.values), max(self.values)
        else:
            vmin, vmax = 0.0, 1.0
        self.lower_spin.setRange(vmin, vmax)
        self.upper_spin.setRange(vmin, vmax)
        self.lower_spin.setDecimals(3)
        self.upper_spin.setDecimals(3)
        self.lower_spin.setValue(vmin)
        self.upper_spin.setValue(vmax)

        self.canvas = HistogramCanvas(self, width=4, height=2, dpi=100)
        # Ensure the plot area grows with the dialog to keep labels visible
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.canvas, stretch=1)

        self.lower_spin.valueChanged.connect(self._update_plot)
        self.upper_spin.valueChanged.connect(self._update_plot)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

        self._update_plot()

    def _update_plot(self) -> None:
        self.canvas.plot_values(
            self.values, self.lower_spin.value(), self.upper_spin.value()
        )

    # -----------------------------
    # Quantile tails support
    # -----------------------------
    def _on_quantile_changed(self, _=None) -> None:
        """Auto-set thresholds from a manual tail percentage.

        Behavior:
        - 0%: do nothing (manual mode).
        - (0, 50%): set thresholds to q and 1-q quantiles.
        - 50%: set both thresholds to the median (0.5 quantile).
        """
        if not self.values:
            return
        pct = float(self.quantile_spin.value())
        # 0% keeps manual thresholds
        if pct <= 0.0:
            return
        # Clamp to [0, 50] for safety (spin box already enforces this)
        pct = min(50.0, max(0.0, pct))
        svals = sorted(float(v) for v in self.values)
        vmin = self.lower_spin.minimum()
        vmax = self.upper_spin.maximum()
        if pct >= 50.0:
            # Both thresholds set to the median
            med = self._percentile(svals, 0.5)
            med = max(vmin, min(vmax, float(med)))
            self.lower_spin.setValue(med)
            self.upper_spin.setValue(med)
            self._update_plot()
            return
        q = pct / 100.0
        low = self._percentile(svals, q)
        high = self._percentile(svals, 1.0 - q)
        low = max(vmin, min(vmax, float(low)))
        high = max(vmin, min(vmax, float(high)))
        if low > high:
            low = high
        self.lower_spin.setValue(low)
        self.upper_spin.setValue(high)
        self._update_plot()

    def _percentile(self, sorted_vals: list[float], q: float) -> float:
        """Return the q-th quantile using linear interpolation (q in [0,1])."""
        n = len(sorted_vals)
        if n == 0:
            return 0.0
        if q <= 0:
            return float(sorted_vals[0])
        if q >= 1:
            return float(sorted_vals[-1])
        pos = q * (n - 1)
        lo = int(pos)
        hi = min(lo + 1, n - 1)
        frac = pos - lo
        return float((1.0 - frac) * sorted_vals[lo] + frac * sorted_vals[hi])

    @property
    def lower_threshold(self) -> float:
        return self.lower_spin.value()

    @property
    def upper_threshold(self) -> float:
        return self.upper_spin.value()


class PercentContrastDialog(QDialog):
    """Dialog to choose minimum percent contrast and inspect threshold-vs-pair sweep."""

    def __init__(
        self,
        sweep_thresholds: Sequence[float],
        sweep_counts: Sequence[int],
        *,
        count_fn: Callable[[float], int] | None = None,
        default_threshold: float | None = None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Percent Contrast Threshold")
        self.setSizeGripEnabled(True)
        self._sweep_thresholds = [float(x) for x in sweep_thresholds]
        self._sweep_counts = [int(x) for x in sweep_counts]
        self._count_fn = count_fn

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        expl = QLabel(
            "Select the minimum % difference for a pair:\n"
            "upper must be at least (1 + %/100) times lower."
        )
        expl.setWordWrap(True)
        layout.addWidget(expl)

        row = QHBoxLayout()
        row.addWidget(QLabel("Minimum % difference:"))
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setDecimals(1)
        self.threshold_spin.setSingleStep(1.0)
        self.threshold_spin.setSuffix("%")
        vmax = max(self._sweep_thresholds) if self._sweep_thresholds else 100.0
        self.threshold_spin.setRange(0.0, max(1.0, float(vmax)))
        if default_threshold is not None:
            self.threshold_spin.setValue(max(0.0, min(float(vmax), float(default_threshold))))
        row.addWidget(self.threshold_spin)
        row.addStretch()
        layout.addLayout(row)

        self.pair_count_lbl = QLabel("")
        layout.addWidget(self.pair_count_lbl)

        self.canvas = HistogramCanvas(self, width=4, height=2, dpi=100)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.canvas, stretch=1)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

        self.threshold_spin.valueChanged.connect(self._update_preview)
        self._update_preview()

    def _nearest_count(self, thr: float) -> int:
        if not self._sweep_thresholds:
            return 0
        best_i = 0
        best_d = abs(self._sweep_thresholds[0] - thr)
        for i in range(1, len(self._sweep_thresholds)):
            d = abs(self._sweep_thresholds[i] - thr)
            if d < best_d:
                best_i = i
                best_d = d
        return int(self._sweep_counts[best_i]) if self._sweep_counts else 0

    def _update_preview(self, _=None) -> None:
        thr = float(self.threshold_spin.value())
        if self._count_fn is not None:
            try:
                count = int(self._count_fn(thr))
            except Exception:
                count = self._nearest_count(thr)
        else:
            count = self._nearest_count(thr)
        self.pair_count_lbl.setText(f"Selected pairs at current threshold: {count}")
        self.canvas.plot_threshold_curve(self._sweep_thresholds, self._sweep_counts, thr)

    @property
    def min_pct_diff(self) -> float:
        return float(self.threshold_spin.value())


class OutgroupDialog(QDialog):
    """Dialog to pick an outgroup species or use ancestral reconstruction."""

    def __init__(self, species: Sequence[str], parent=None, default_selected: str | None = None, show_two_pair_option: bool = False, default_agreement_pct: float | None = 100.0, species_groups_file: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Select Outgroup Species")
        self._species_groups_file = species_groups_file
        layout = QVBoxLayout(self)
        
        # Radio button group for outgroup method
        method_label = QLabel("<b>Outgroup Method:</b>")
        layout.addWidget(method_label)
        
        self._method_group = QButtonGroup(self)
        self._use_species_radio = QRadioButton("Use single outgroup species")
        self._use_ancestral_radio = QRadioButton("Use parsimony ancestral reconstruction")
        self._use_species_radio.setChecked(True)
        
        self._method_group.addButton(self._use_species_radio, 0)
        self._method_group.addButton(self._use_ancestral_radio, 1)
        
        layout.addWidget(self._use_species_radio)
        layout.addWidget(self._use_ancestral_radio)
        
        # Species selection (for single species method)
        species_label = QLabel("Choose outgroup species:")
        layout.addWidget(species_label)
        self._species_label = species_label
        
        self._combo = QComboBox()
        self._combo.addItems(list(species))
        # Preselect default if provided and present in the list
        if default_selected:
            idx = self._combo.findText(default_selected)
            if idx >= 0:
                self._combo.setCurrentIndex(idx)
        layout.addWidget(self._combo)
        
        # Tree file selection (for ancestral method)
        tree_label = QLabel("Tree file (Newick/NEXUS):")
        layout.addWidget(tree_label)
        self._tree_label = tree_label
        
        tree_row = QHBoxLayout()
        self._tree_path = QLineEdit()
        self._tree_path.setPlaceholderText("Select tree file for ancestral reconstruction...")
        self._tree_path.setReadOnly(True)
        tree_row.addWidget(self._tree_path)
        
        self._tree_browse_btn = QPushButton("Browse...")
        self._tree_browse_btn.clicked.connect(self._browse_tree)
        tree_row.addWidget(self._tree_browse_btn)
        
        layout.addLayout(tree_row)
        
        # Validation label for tree
        self._tree_validation_label = QLabel()
        self._tree_validation_label.setWordWrap(True)
        self._tree_validation_label.setStyleSheet("QLabel { color: gray; font-style: italic; }")
        layout.addWidget(self._tree_validation_label)
        
        # Connect radio buttons to update UI
        self._use_species_radio.toggled.connect(self._update_method_ui)
        
        # Initial UI state
        self._update_method_ui()
        # Option: use 2x2 combos derived from species groups (Fast Scan only)
        self._two_pair_check = None
        if show_two_pair_option:
            self._two_pair_check = QCheckBox("Use all two-pair combos")
            self._two_pair_check.setToolTip(
                "Interpret the species groups file by running every combination of two pairs (2x2)\n"
                "instead of the standard NxN across all pairs. For pairs with multiple species\n"
                "options on a line, generate variants for each option."
            )
            layout.addWidget(self._two_pair_check)
        # Minimum outgroup-control agreement percentage (0–100%)
        agree_row = QHBoxLayout()
        agree_row.addWidget(QLabel("Minimum outgroup-control agreement (%):"))
        self._agree_spin = QDoubleSpinBox()
        self._agree_spin.setSuffix("%")
        self._agree_spin.setRange(0.0, 100.0)
        self._agree_spin.setDecimals(0)
        self._agree_spin.setSingleStep(5.0)
        try:
            init_pct = float(default_agreement_pct if default_agreement_pct is not None else 100.0)
        except Exception:
            init_pct = 100.0
        init_pct = max(0.0, min(100.0, init_pct))
        self._agree_spin.setValue(init_pct)
        agree_row.addWidget(self._agree_spin)
        agree_row.addStretch()
        layout.addLayout(agree_row)
        
        # Require unambiguous MRCA checkbox (for ancestral mode)
        self._require_unamb_check = QCheckBox("Require unambiguous ancestral residues")
        self._require_unamb_check.setToolTip(
            "When checked, sites with ambiguous ancestral states (multiple equally parsimonious residues)\n"
            "are excluded from CCS detection. By default (unchecked), ambiguous sites are evaluated\n"
            "using the set of possible ancestral states."
        )
        self._require_unamb_check.setChecked(False)
        layout.addWidget(self._require_unamb_check)
        
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

        # Ensure the dialog is initially wide enough so the window title isn't truncated
        try:
            min_w = 420
            self.setMinimumWidth(min_w)
            # Use sizeHint to get a sensible initial height; ensure width >= min_w
            hint = self.sizeHint()
            self.resize(max(hint.width(), min_w), hint.height())
        except Exception:
            pass

    @property
    def selected(self) -> str:
        return self._combo.currentText()

    @property
    def use_two_pair_combos(self) -> bool:
        try:
            return bool(self._two_pair_check and self._two_pair_check.isChecked())
        except Exception:
            return False

    @property
    def min_out_ctrl_agreement(self) -> float:
        """Return minimum outgroup-control agreement as a fraction [0,1]."""
        try:
            return float(self._agree_spin.value()) / 100.0
        except Exception:
            return 1.0
    
    @property
    def use_ancestral_reconstruction(self) -> bool:
        """Return True if ancestral reconstruction is selected."""
        return self._use_ancestral_radio.isChecked()
    
    @property
    def tree_file(self) -> str:
        """Return path to tree file (empty if not using ancestral reconstruction)."""
        return self._tree_path.text().strip()
    
    def _browse_tree(self) -> None:
        """Open file dialog to select tree file."""
        import os
        start_dir = os.getcwd()
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Tree File",
            start_dir,
            "Tree Files (*.nwk *.newick *.tree *.tre *.treefile *.contree *.nh *.nhx *.dnd *.nexus *.nex *.txt);;All Files (*)"
        )
        if path:
            self._tree_path.setText(path)
            self._validate_tree()
    
    def _validate_tree(self) -> None:
        """Validate selected tree against species groups."""
        tree_path = self._tree_path.text().strip()
        if not tree_path:
            self._tree_validation_label.setText("")
            self._tree_validation_label.setStyleSheet("QLabel { color: gray; font-style: italic; }")
            return
        
        if not self._species_groups_file:
            self._tree_validation_label.setText("⚠ Cannot validate: no species groups file loaded")
            self._tree_validation_label.setStyleSheet("QLabel { color: orange; font-style: italic; }")
            return
        
        try:
            from gui.core.ancestral_reconstruction import validate_tree_for_fast_scan
            valid, message = validate_tree_for_fast_scan(tree_path, self._species_groups_file)
            if valid:
                self._tree_validation_label.setText(f"✓ {message}")
                self._tree_validation_label.setStyleSheet("QLabel { color: green; }")
            else:
                self._tree_validation_label.setText(f"✗ {message}")
                self._tree_validation_label.setStyleSheet("QLabel { color: red; }")
        except Exception as e:
            self._tree_validation_label.setText(f"✗ Validation failed: {e}")
            self._tree_validation_label.setStyleSheet("QLabel { color: red; }")
    
    def _update_method_ui(self) -> None:
        """Update UI based on selected outgroup method."""
        use_species = self._use_species_radio.isChecked()
        
        # Show/hide species controls
        self._species_label.setVisible(use_species)
        self._combo.setVisible(use_species)
        
        # Show/hide tree controls
        self._tree_label.setVisible(not use_species)
        self._tree_path.setVisible(not use_species)
        self._tree_browse_btn.setVisible(not use_species)
        self._tree_validation_label.setVisible(not use_species)
        # Show 'require unambiguous' only for ancestral mode
        try:
            self._require_unamb_check.setVisible(not use_species)
        except Exception:
            pass
        
        # Validate tree if switching to ancestral mode
        if not use_species:
            self._validate_tree()
    
    def require_unambiguous_mrca(self) -> bool:
        """Return whether user requires unambiguous ancestral residues."""
        return self._require_unamb_check.isChecked()
