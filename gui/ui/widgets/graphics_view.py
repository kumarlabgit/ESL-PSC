from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import QGraphicsView, QGraphicsTextItem
from PySide6.QtGui import QPainter, QPen, QColor, QLinearGradient, QBrush, QPalette
from PySide6.QtCore import Qt

from gui.ui.widgets.color_utils import viridis_stops_rgb


class ZoomableGraphicsView(QGraphicsView):
    """Graphics view that supports wheel-based zooming with limits and draws a
    small Viridis color legend when the parent provides continuous phenotype data.

    This class avoids importing TreeViewer directly to prevent circular imports.
    It uses duck typing to interact with its parent:
    - context menu routing if parent has `_show_label_menu` / `_show_pair_menu`
    - draws legend if parent has `_continuous_pheno` set to True
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._zoom_level = 0
        self._min_zoom = -5
        self._max_zoom = 10
        # Prevent smearing of foreground overlays (like the legend) during
        # interactive panning/zooming by forcing full viewport repaints.
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

    # ------------------------------------------------------------------
    def wheelEvent(self, event):
        if event.angleDelta().y() > 0 and self._zoom_level < self._max_zoom:
            factor = 1.25
            self._zoom_level += 1
            self.scale(factor, factor)
        elif event.angleDelta().y() < 0 and self._zoom_level > self._min_zoom:
            factor = 0.8
            self._zoom_level -= 1
            self.scale(factor, factor)

    # ------------------------------------------------------------------
    def resetTransform(self):
        super().resetTransform()
        self._zoom_level = 0

    # ------------------------------------------------------------------
    def contextMenuEvent(self, event):
        parent = self.parentWidget()
        # Route to parent's label/pair menus if available; otherwise fallback.
        if parent is None:
            super().contextMenuEvent(event)
            return
        items = self.items(event.pos())
        label_item: Optional[QGraphicsTextItem] = None
        pair_item: Optional[QGraphicsTextItem] = None
        for it in items:
            if isinstance(it, QGraphicsTextItem) and hasattr(it, "species_name"):
                label_item = it
                break
            if isinstance(it, QGraphicsTextItem) and hasattr(it, "pair_index"):
                pair_item = it
        if label_item is not None and hasattr(parent, "_show_label_menu"):
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            label_item.setCursor(Qt.CursorShape.OpenHandCursor)
            parent._show_label_menu(label_item, self.mapToGlobal(event.pos()))
        elif pair_item is not None and hasattr(parent, "_show_pair_menu"):
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            pair_item.setCursor(Qt.CursorShape.OpenHandCursor)
            parent._show_pair_menu(pair_item, self.mapToGlobal(event.pos()))
        else:
            super().contextMenuEvent(event)

    # ------------------------------------------------------------------
    def drawForeground(self, painter: QPainter, rect):
        """Draw a small viridis color bar legend in the view's upper-left.

        Anchored in view coordinates so it remains unobtrusive regardless of
        scene zoom/pan. Only shown when parent indicates continuous phenotypes.
        """
        super().drawForeground(painter, rect)
        parent = self.parentWidget()
        if not (parent is not None and getattr(parent, "_continuous_pheno", False)):
            return

        # Switch to view (pixel) coordinates
        painter.save()
        painter.resetTransform()

        margin = 8
        bar_w = 120
        bar_h = 10

        pal = self.palette()
        text_color = pal.color(QPalette.ColorRole.WindowText)
        painter.setPen(text_color)

        # Title above the bar
        title = "trait values"
        fm = painter.fontMetrics()
        title_x = margin
        title_y = margin + fm.ascent()
        painter.drawText(title_x, title_y, title)

        # Color bar rectangle just below the title
        bar_x = margin
        bar_y = margin + fm.height() + 2

        grad = QLinearGradient(bar_x, bar_y, bar_x + bar_w, bar_y)
        for t, (r, g, b) in viridis_stops_rgb():
            grad.setColorAt(t, QColor(r, g, b))

        painter.setBrush(QBrush(grad))
        painter.setPen(QPen(text_color, 1))
        painter.drawRect(bar_x, bar_y, bar_w, bar_h)

        # End labels under the bar
        low_txt = "low"
        high_txt = "high"
        labels_y = bar_y + bar_h + 2 + fm.ascent()
        painter.setPen(text_color)
        painter.drawText(bar_x, labels_y, low_txt)
        hi_w = fm.horizontalAdvance(high_txt)
        painter.drawText(bar_x + bar_w - hi_w, labels_y, high_txt)

        painter.restore()
