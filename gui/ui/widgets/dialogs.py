from __future__ import annotations

from typing import Sequence

from PySide6.QtWidgets import (
    QMessageBox,
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QDialogButtonBox,
    QLabel,
    QDoubleSpinBox,
)

from gui.ui.widgets.histogram_canvas import HistogramCanvas


class AutoSelectOptionsDialog(QMessageBox):
    """Simple dialog to choose auto-select behavior.

    Presents options in this order: Longest Sequence → Max trait contrast → Random → Default → Cancel.
    Stores the chosen option as `self.choice` in {"default","longest","random","contrast",None}.
    """

    def __init__(self, allow_longest: bool, allow_contrast: bool, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Auto Select Options")
        self.setText("Choose how to resolve equally valid siblings")

        # Buttons in a fixed order using ActionRole to preserve ordering cross-platform
        self.longest_btn = self.addButton(
            "Longest Sequence", QMessageBox.ButtonRole.ActionRole
        )
        if not allow_longest:
            self.longest_btn.setEnabled(False)
            self.longest_btn.setToolTip(
                "Set an alignments directory in the input page, to enable this option"
            )

        self.contrast_btn = self.addButton(
            "Max trait contrast", QMessageBox.ButtonRole.ActionRole
        )
        if not allow_contrast:
            self.contrast_btn.setEnabled(False)
            self.contrast_btn.setToolTip(
                "Available only for continuous phenotypes"
            )

        self.random_btn = self.addButton(
            "Random", QMessageBox.ButtonRole.ActionRole
        )

        self.default_btn = self.addButton(
            "Default", QMessageBox.ButtonRole.ActionRole
        )
        # Make the Default button appear as the default (blue) accept button.
        self.setDefaultButton(self.default_btn)

        # Cancel button
        self.cancel_btn = self.addButton(
            "Cancel", QMessageBox.ButtonRole.RejectRole
        )

        self.choice: str | None = None

    def exec(self) -> int:
        result = super().exec()
        clicked = self.clickedButton()
        if clicked == self.default_btn:
            self.choice = "default"
        elif clicked == self.longest_btn:
            self.choice = "longest"
        elif clicked == self.contrast_btn:
            self.choice = "contrast"
        elif clicked == self.random_btn:
            self.choice = "random"
        else:
            self.choice = None
        return result

    def reject(self) -> None:
        """Handle programmatic rejection (e.g., Esc)."""
        self.choice = None
        super().reject()

    def closeEvent(self, event):
        """Treat window close as cancel."""
        self.choice = None
        super().closeEvent(event)


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

    @property
    def lower_threshold(self) -> float:
        return self.lower_spin.value()

    @property
    def upper_threshold(self) -> float:
        return self.upper_spin.value()
