from PyQt6.QtWidgets import QWidget
from PyQt6.QtGui import QPainter, QColor, QPen


class ProteinMapWidget(QWidget):
    """Simple schematic showing selected site positions on a protein."""

    def __init__(self, length: int, sites: list[int], parent=None):
        super().__init__(parent)
        self.length = max(length, 1)
        self.sites = sites
        self.setMinimumHeight(20)
        self.setMaximumHeight(20)

    def paintEvent(self, event):
        painter = QPainter(self)
        bar_height = 6
        y = (self.height() - bar_height) // 2
        painter.fillRect(0, y, self.width(), bar_height, QColor(220, 220, 220))
        pen = QPen(QColor("green"))
        for site in self.sites:
            x = int(site / self.length * self.width())
            painter.setPen(pen)
            painter.drawLine(x, y, x, y + bar_height)
        painter.end()
