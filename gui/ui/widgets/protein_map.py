from __future__ import annotations

from PySide6.QtWidgets import QWidget
from PySide6.QtGui import QPainter, QColor, QPen
import math


class ProteinMapWidget(QWidget):
    """Simple schematic showing selected site positions on a protein."""

    def __init__(self, length: int, sites: list[int], scores: list[float] | None = None, parent=None, score_scale_max: float | None = None):
        super().__init__(parent)
        self.length = max(length, 1)
        self.sites = sites
        self.scores = scores or []
        # Optional global scale for opacity across multiple proteins:
        # when provided, 0.0 -> fully transparent and score_scale_max -> fully opaque.
        self.score_scale_max = score_scale_max
        self.setMinimumHeight(30)  # Increased minimum height for better visibility
        # Ensure the widget is wide enough to be visible when embedded in tables
        self.setMinimumWidth(500)
        self.setMaximumHeight(50)  # Increased maximum height to accommodate extended lines

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Draw protein bar
        bar_height = 10  # Increased bar height
        y = (self.height() - bar_height) // 2
        
        # Use a color that works well in both light and dark modes
        bg_color = QColor(150, 150, 150) if self.palette().window().color().lightness() < 128 else QColor(220, 220, 220)
        painter.fillRect(0, y, self.width(), bar_height, bg_color)
        
        # Configure pen for site lines
        # Pre-compute score scaling if scores are provided.
        # If a global max was supplied, scale transparency using a logarithmic mapping
        # relative to [0, score_scale_max]: rel = log1p(val) / log1p(score_scale_max).
        min_score, max_score = None, None
        if self.scores:
            if self.score_scale_max is not None and self.score_scale_max > 0:
                min_score = 0.0
                max_score = float(self.score_scale_max)
            else:
                min_score = min(self.scores)
                max_score = max(self.scores)
                if max_score == min_score:
                    max_score += 1e-9  # avoid division by zero

        # Draw site lines extending above and below the protein bar
        extension = 3  # Pixels to extend lines above and below
        for idx, site in enumerate(self.sites):
            if self.scores and idx < len(self.scores):
                val = self.scores[idx]
                if self.score_scale_max is not None and self.score_scale_max > 0:
                    # Logarithmic scaling across the global max.
                    # Map val<=0 to 0 to avoid log domain issues.
                    val = max(0.0, float(val))
                    denom = math.log1p(max_score) if max_score is not None else 0.0
                    if denom > 0.0:
                        rel = math.log1p(val) / denom
                    else:
                        rel = 0.0
                else:
                    # Legacy per-protein linear scaling
                    rel = (val - min_score) / (max_score - min_score)
                if rel < 0.0:
                    rel = 0.0
                elif rel > 1.0:
                    rel = 1.0
            else:
                rel = 1.0  # default full visibility when scores missing

            # Base color (same hue as original #4CAF50)
            base_color = QColor("#4CAF50")
            # Blend toward background color as scores decrease
            r = int(base_color.red() * rel + bg_color.red() * (1 - rel))
            g = int(base_color.green() * rel + bg_color.green() * (1 - rel))
            b = int(base_color.blue() * rel + bg_color.blue() * (1 - rel))
            alpha = int(255 * rel)
            pen_color = QColor(r, g, b, alpha)
            pen = QPen(pen_color)
            pen.setWidth(2)
            painter.setPen(pen)

            x = int(site / self.length * self.width())
            painter.drawLine(x, y - extension, x, y + bar_height + extension)

        painter.end()
