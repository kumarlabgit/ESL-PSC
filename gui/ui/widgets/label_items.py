from __future__ import annotations

from typing import Optional, Dict

from PySide6.QtWidgets import QGraphicsTextItem
from PySide6.QtCore import Qt


class HoverLabelItem(QGraphicsTextItem):
    """Text item that changes cursor on hover and supports quick phenotype assign.

    This implementation avoids importing TreeViewer to prevent circular imports.
    It uses duck-typing against its parent (the widget that owns the view):
    - expects optional attributes: _assign_mode, _assign_cursors, _assign_pheno
    - reads a dynamic attribute 'species_name' set by the caller
    """

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.OpenHandCursor)

    # ------------------------------------------------------------------
    def _parent_widget(self):
        views = self.scene().views() if self.scene() else []
        return views[0].parentWidget() if views else None

    # ------------------------------------------------------------------
    def hoverEnterEvent(self, event):
        parent = self._parent_widget()
        mode = getattr(parent, "_assign_mode", None)
        cursors: Dict[int, object] = getattr(parent, "_assign_cursors", {}) if parent else {}
        if mode in (1, -1):
            cur = cursors.get(mode)
            if cur is not None:
                self.setCursor(cur)
            else:
                self.setCursor(Qt.CursorShape.ArrowCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        super().hoverEnterEvent(event)

    # ------------------------------------------------------------------
    def hoverLeaveEvent(self, event):
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        super().hoverLeaveEvent(event)

    # ------------------------------------------------------------------
    def mousePressEvent(self, event):
        parent = self._parent_widget()
        if (
            event.button() == Qt.MouseButton.LeftButton
            and parent is not None
            and hasattr(self, "species_name")
        ):
            mode = getattr(parent, "_assign_mode", None)
            assign_fn = getattr(parent, "_assign_pheno", None)
            if mode in (1, -1) and callable(assign_fn):
                assign_fn(self.species_name, mode)
                event.accept()
                return
        super().mousePressEvent(event)
