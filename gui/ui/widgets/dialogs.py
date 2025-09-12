from __future__ import annotations

from typing import Sequence

from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QDialogButtonBox,
    QLabel,
    QDoubleSpinBox,
    QPushButton,
    QComboBox,
)

from gui.ui.widgets.histogram_canvas import HistogramCanvas


class AutoSelectOptionsDialog(QDialog):
    """Dialog to choose auto-select behavior with vertically stacked buttons.

    Presents options in this order: Longest Sequence → Shortest distance → Max trait contrast → Composite best → Random → Default → Cancel.
    Stores the chosen option as `self.choice` in {"default","longest","shortest","contrast","composite","random",None}.
    """

    def __init__(self, allow_longest: bool, allow_contrast: bool, allow_composite: bool = True, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Auto Select Options")
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Choose how to resolve equally valid siblings"))

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

    def _pick(self, choice: str) -> None:
        self.choice = choice
        self.accept()


class PhenoThresholdDialog(QDialog):
    """Dialog for choosing thresholds for continuous phenotypes."""

    def __init__(self, values: Sequence[float], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phenotype Thresholds")
        layout = QVBoxLayout(self)

        form = QHBoxLayout()
        form.addWidget(QLabel("Lower threshold:"))
        self.lower_spin = QDoubleSpinBox()
        form.addWidget(self.lower_spin)
        form.addWidget(QLabel("Upper threshold:"))
        self.upper_spin = QDoubleSpinBox()
        form.addWidget(self.upper_spin)
        layout.addLayout(form)

        # Quantile tail auto-set controls
        qrow = QHBoxLayout()
        qlbl = QLabel("Auto-set by tails:")
        qlbl.setToolTip(
            "Choose a tail percentage to automatically set thresholds to the lower and upper quantiles\n"
            "(e.g., 10% sets thresholds at the 10th and 90th percentiles)."
        )
        self.quantile_combo = QComboBox()
        # First entry keeps manual control; subsequent entries auto-apply when selected
        self.quantile_combo.addItems(["Manual", "5%", "10%", "15%", "20%", "25%", "30%"])
        self.quantile_combo.setToolTip(
            "Automatically set thresholds using symmetric quantile tails."
        )
        self.quantile_combo.currentIndexChanged.connect(self._on_quantile_changed)
        qrow.addWidget(qlbl)
        qrow.addWidget(self.quantile_combo)
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
        layout.addWidget(self.canvas)

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
    def _on_quantile_changed(self, idx: int) -> None:
        """When a quantile tail is chosen, auto-set lower/upper to matching percentiles."""
        if not self.values:
            return
        text = self.quantile_combo.currentText().strip()
        if not text or text.lower().startswith("manual"):
            return  # keep manual values
        # Parse like "10%"
        try:
            pct = float(text.rstrip("%"))
        except Exception:
            return
        if pct <= 0 or pct >= 50:
            # Tails must be in (0, 50); ignore otherwise
            return
        q = pct / 100.0
        svals = sorted(float(v) for v in self.values)
        low = self._percentile(svals, q)
        high = self._percentile(svals, 1.0 - q)
        # Apply with clamping to spin ranges
        vmin = self.lower_spin.minimum()
        vmax = self.upper_spin.maximum()
        low = max(vmin, min(vmax, float(low)))
        high = max(vmin, min(vmax, float(high)))
        if low > high:
            # Degenerate distribution; force equality
            low = high
        # Set without re-triggering quantile change logic
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


class OutgroupDialog(QDialog):
    """Dialog to pick an outgroup species."""

    def __init__(self, species: Sequence[str], parent=None, default_selected: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Select Outgroup Species")
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Choose outgroup species:"))
        self._combo = QComboBox()
        self._combo.addItems(list(species))
        # Preselect default if provided and present in the list
        if default_selected:
            idx = self._combo.findText(default_selected)
            if idx >= 0:
                self._combo.setCurrentIndex(idx)
        layout.addWidget(self._combo)
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    @property
    def selected(self) -> str:
        return self._combo.currentText()
