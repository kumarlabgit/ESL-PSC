from __future__ import annotations

import math
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

    _ANGLE_UNITS_PER_STEP = 120
    _PIXELS_PER_TRACKPAD_STEP = 900.0
    _MAX_TRACKPAD_STEPS_PER_EVENT = 0.12
    _MAX_WHEEL_STEPS_PER_EVENT = 3

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._zoom_level = 0.0
        self._min_zoom = -5
        self._max_zoom = 10
        self._angle_zoom_remainder = 0
        # Prevent smearing of foreground overlays (like the legend) during
        # interactive panning/zooming by forcing full viewport repaints.
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

    # ------------------------------------------------------------------
    def wheelEvent(self, event):
        pixel_delta = event.pixelDelta().y()
        if pixel_delta:
            # Trackpads report high-frequency pixel deltas. Keep each event
            # small so macOS two-finger gestures do not jump through the tree.
            steps = max(
                -self._MAX_TRACKPAD_STEPS_PER_EVENT,
                min(self._MAX_TRACKPAD_STEPS_PER_EVENT, pixel_delta / self._PIXELS_PER_TRACKPAD_STEP),
            )
            self._apply_zoom_steps(steps)
            event.accept()
            return

        angle_delta = event.angleDelta().y()
        if not angle_delta:
            super().wheelEvent(event)
            return

        self._angle_zoom_remainder += angle_delta
        steps = int(self._angle_zoom_remainder / self._ANGLE_UNITS_PER_STEP)
        if steps == 0:
            event.accept()
            return

        self._angle_zoom_remainder -= steps * self._ANGLE_UNITS_PER_STEP
        self._apply_zoom_steps(max(-self._MAX_WHEEL_STEPS_PER_EVENT, min(self._MAX_WHEEL_STEPS_PER_EVENT, steps)))
        event.accept()

    # ------------------------------------------------------------------
    def _apply_zoom_steps(self, steps: float) -> None:
        if steps == 0:
            return
        target = max(self._min_zoom, min(self._max_zoom, self._zoom_level + steps))
        actual_steps = target - self._zoom_level
        if actual_steps == 0:
            return
        factor = math.pow(1.15, actual_steps)
        self._zoom_level = target
        self.scale(factor, factor)

    # ------------------------------------------------------------------
    def resetTransform(self):
        super().resetTransform()
        self._zoom_level = 0.0
        self._angle_zoom_remainder = 0

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
