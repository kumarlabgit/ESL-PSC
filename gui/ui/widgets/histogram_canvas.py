from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class HistogramCanvas(FigureCanvasQTAgg):
    """Matplotlib canvas for displaying a simple histogram."""

    def __init__(self, parent=None, width: float = 4, height: float = 2, dpi: int = 100) -> None:
        fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(fig)
        self.setParent(parent)
        self.ax = fig.add_subplot(111)

    def plot_scores(self, scores: list[float], threshold: float) -> None:
        """Plot a histogram of scores and draw a threshold line."""
        self.ax.clear()
        if scores:
            self.ax.hist(scores, bins=20, color="lightblue", edgecolor="black")
            self.ax.axvline(threshold, color="red", linestyle="--")
        self.ax.set_xlabel("Score")
        self.draw()
