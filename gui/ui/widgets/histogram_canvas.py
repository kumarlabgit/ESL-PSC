from __future__ import annotations

from typing import Sequence

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import math
import numpy as np
from PySide6.QtWidgets import QSizePolicy


class HistogramCanvas(FigureCanvasQTAgg):
    """Matplotlib histogram canvas usable in PySide6 widgets.

    Uses constrained_layout so that axes, tick labels, and axis labels are
    automatically laid out to fit within the canvas bounds at any size.
    """

    def __init__(self, parent=None, width: float = 5, height: float = 2, dpi: int = 100) -> None:
        fig = Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        super().__init__(fig)
        self.setParent(parent)
        # Make the canvas expand to fill available space in layouts
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
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

    def plot_values(
        self, values: Sequence[float], lower: float, upper: float
    ) -> None:
        """Plot a histogram of phenotype values with threshold lines."""
        self.axes.clear()

        if not values:
            self.draw()
            return

        # Decide whether to use a log-scaled x-axis.
        # Criteria: all values are strictly positive AND span >= 2 orders of magnitude.
        vals = np.asarray(list(values), dtype=float)
        pos_mask = vals > 0
        use_log = False
        if np.all(pos_mask):
            vmin = float(np.min(vals))
            vmax = float(np.max(vals))
            if vmin > 0 and vmax > 0:
                orders = (math.log10(vmax) - math.log10(vmin)) if vmax > vmin else 0.0
                use_log = orders >= 2.0

        if use_log:
            # Log-spaced bins over the positive range.
            vmin = float(np.min(vals))
            vmax = float(np.max(vals))
            # Choose a sensible number of bins based on dynamic range (cap at 60)
            # e.g., ~15 bins per decade, minimum 20 bins overall
            orders = max(1.0, math.log10(vmax) - math.log10(vmin))
            num_bins = int(min(60, max(20, 15 * orders)))
            bins = np.logspace(math.log10(vmin), math.log10(vmax), num_bins)
            self.axes.hist(
                vals,
                bins=bins,
                color="lightblue",
                edgecolor="black",
            )
            self.axes.set_xscale("log")
            # Draw threshold lines only if thresholds are positive in log mode
            if lower > 0:
                self.axes.axvline(
                    x=lower,
                    color="red",
                    linestyle="--",
                    label=f"Lower = {lower:.3g}",
                )
            if upper > 0:
                self.axes.axvline(
                    x=upper,
                    color="green",
                    linestyle="--",
                    label=f"Upper = {upper:.3g}",
                )
            self.axes.set(
                xlabel="Phenotype Value (log scale)",
                ylabel="Count of Species",
                title="Phenotype Value Distribution",
            )
        else:
            # Fallback to linear scale
            self.axes.hist(
                vals,
                bins="auto",
                color="lightblue",
                edgecolor="black",
            )
            self.axes.axvline(
                x=lower,
                color="red",
                linestyle="--",
                label=f"Lower = {lower:.3g}",
            )
            self.axes.axvline(
                x=upper,
                color="green",
                linestyle="--",
                label=f"Upper = {upper:.3g}",
            )
            self.axes.set(
                xlabel="Phenotype Value",
                ylabel="Count of Species",
                title="Phenotype Value Distribution",
            )
        self.axes.legend()
        # With constrained_layout enabled, drawing triggers a layout pass that
        # allocates enough space for tick/axis labels to avoid clipping.
        self.draw()

    def plot_threshold_curve(
        self,
        thresholds: Sequence[float],
        pair_counts: Sequence[int],
        selected_threshold: float | None = None,
    ) -> None:
        """Plot threshold vs selected-pair count for percent-contrast tuning."""
        self.axes.clear()
        xs = [float(x) for x in thresholds]
        ys = [int(y) for y in pair_counts]
        if not xs or not ys or len(xs) != len(ys):
            self.draw()
            return
        self.axes.plot(xs, ys, marker="o", markersize=3, color="#1f77b4")
        if selected_threshold is not None:
            self.axes.axvline(
                x=float(selected_threshold),
                color="red",
                linestyle="--",
                label=f"Selected = {float(selected_threshold):.1f}%",
            )
            self.axes.legend()
        self.axes.set(
            xlabel="Minimum % difference (upper vs lower)",
            ylabel="Selected pairs",
            title="Threshold Sweep",
        )
        self.draw()

    # Ensure that on resize, a new draw is scheduled so constrained_layout can
    # adapt the subplot geometry to the new size and keep labels visible.
    def resizeEvent(self, event) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        self.draw_idle()
