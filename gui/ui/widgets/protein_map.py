from __future__ import annotations

from PyQt6.QtWidgets import QWidget
from PyQt6.QtGui import QPainter, QColor, QPen


class ProteinMapWidget(QWidget):
    """Simple schematic showing selected site positions on a protein."""

    def __init__(self, length: int, sites: list[int], scores: list[float] | None = None, parent=None):
        super().__init__(parent)
        self.length = max(length, 1)
        self.sites = sites
        self.scores = scores or []
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
        # Pre-compute score scaling if scores are provided
        min_score, max_score = None, None
        if self.scores:
            min_score = min(self.scores)
            max_score = max(self.scores)
            if max_score == min_score:
                max_score += 1e-9  # avoid division by zero

        # Draw site lines extending above and below the protein bar
        extension = 3  # Pixels to extend lines above and below
        for idx, site in enumerate(self.sites):
            # Determine brightness factor based on relative score
            if self.scores and idx < len(self.scores):
                rel = (self.scores[idx] - min_score) / (max_score - min_score)
                # Keep brightness in [0.2, 1.0] so extremely low scores are still visible faintly
                brightness = 0.2 + 0.8 * rel
            else:
                brightness = 1.0  # default full brightness when scores missing

            # Base color (same hue as original #4CAF50)
            base_color = QColor("#4CAF50")
            r = int(base_color.red() * brightness)
            g = int(base_color.green() * brightness)
            b = int(base_color.blue() * brightness)
            pen = QPen(QColor(r, g, b))
            pen.setWidth(2)
            painter.setPen(pen)

            x = int(site / self.length * self.width())
            painter.drawLine(x, y - extension, x, y + bar_height + extension)
            
        painter.end()
