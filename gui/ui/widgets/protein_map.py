from PyQt6.QtWidgets import QWidget
from PyQt6.QtGui import QPainter, QColor, QPen


class ProteinMapWidget(QWidget):
    """Simple schematic showing selected site positions on a protein."""

    def __init__(self, length: int, sites: list[int], parent=None):
        super().__init__(parent)
        self.length = max(length, 1)
        self.sites = sites
        self.setMinimumHeight(30)  # Increased minimum height for better visibility
        # Ensure the widget is wide enough to be visible when embedded in tables
        self.setMinimumWidth(500)
        self.setMaximumHeight(60)  # Increased maximum height to accommodate extended lines

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
        pen = QPen(QColor("#4CAF50"))  # Brighter green for better visibility
        pen.setWidth(2)  # Thicker lines
        painter.setPen(pen)
        
        # Draw site lines extending above and below the protein bar
        extension = 5  # Pixels to extend lines above and below
        for site in self.sites:
            x = int(site / self.length * self.width())
            # Draw line from top extension to bottom extension through the bar
            painter.drawLine(x, y - extension, x, y + bar_height + extension)
            
        painter.end()
