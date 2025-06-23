from __future__ import annotations

from typing import Sequence

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class HistogramCanvas(FigureCanvasQTAgg):
    """Matplotlib histogram canvas usable in PyQt6 widgets."""

    def __init__(self, parent=None, width: float = 5, height: float = 2, dpi: int = 100) -> None:
        fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(fig)
        self.setParent(parent)
        self.axes = fig.add_subplot(111)

    def plot_scores(self, scores: Sequence[float], threshold: float) -> None:
        """Redraw the histogram with given scores and threshold."""
        self.axes.clear()

        if not scores:
            self.draw()
            return

        b_min, b_max = min(scores), max(scores)
        if b_min == b_max:
            b_min -= 1
            b_max += 1

        # Use one bin per integer score value for discrete look
        self.axes.hist(
            scores,
            bins=(b_max - b_min + 1),
            range=(b_min - 0.5, b_max + 0.5),
            color="lightblue",
            edgecolor="black",
        )
        self.axes.axvline(x=threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")
        self.axes.set(
            xlabel="Convergence Score",
            ylabel="Count of Sites",
            title="Distribution of Convergence Scores",
        )
        self.axes.legend()
        self.draw()
