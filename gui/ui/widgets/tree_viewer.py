from __future__ import annotations

"""A simple QGraphicsView-based viewer for phylogenetic trees."""

from typing import Dict
from PyQt6.QtWidgets import (
    QWidget, QGraphicsView, QGraphicsScene, QVBoxLayout, QGraphicsTextItem
)
from PyQt6.QtGui import QPainter, QPen
from PyQt6.QtCore import Qt
from Bio.Phylo.Newick import Tree, Clade


class _ZoomableGraphicsView(QGraphicsView):
    """Graphics view that supports wheel-based zooming."""

    def wheelEvent(self, event):
        factor = 1.25 if event.angleDelta().y() > 0 else 0.8
        self.scale(factor, factor)


class TreeViewer(QWidget):
    """Window displaying a Newick phylogenetic tree."""

    def __init__(self, tree: Tree, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phylogenetic Tree")
        layout = QVBoxLayout(self)

        self.view = _ZoomableGraphicsView()
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.view.setTransformationAnchor(
            QGraphicsView.ViewportAnchor.AnchorUnderMouse
        )

        layout.addWidget(self.view)
        self.scene = QGraphicsScene(self)
        self.view.setScene(self.scene)

        self._draw_tree(tree)

    # ------------------------------------------------------------------
    def _y_positions(self, tree: Tree, step: int = 30) -> Dict[Clade, float]:
        """Assign y positions to all nodes based on tip order."""
        y: Dict[Clade, float] = {}
        for idx, leaf in enumerate(tree.get_terminals()):
            y[leaf] = idx * step

        def set_internal(clade: Clade) -> float:
            if clade.is_terminal():
                return y[clade]
            vals = [set_internal(c) for c in clade.clades]
            y[clade] = sum(vals) / len(vals)
            return y[clade]

        set_internal(tree.root)
        return y

    # ------------------------------------------------------------------
    def _draw_tree(self, tree: Tree) -> None:
        """Render the tree to the graphics scene."""
        depths = tree.depths(unit_branch_lengths=True)
        max_x = max(depths.values()) if depths else 0
        y_pos = self._y_positions(tree)

        pen = QPen(Qt.GlobalColor.black)
        for clade in tree.find_clades(order="preorder"):
            x1, y1 = depths.get(clade, 0), y_pos.get(clade, 0)
            for child in clade.clades:
                x2, y2 = depths.get(child, 0), y_pos.get(child, 0)
                self.scene.addLine(x1, y1, x2, y2, pen)

        # horizontal lines to a common x
        for leaf in tree.get_terminals():
            x_leaf = depths.get(leaf, 0)
            y_leaf = y_pos.get(leaf, 0)
            self.scene.addLine(x_leaf, y_leaf, max_x, y_leaf, pen)
            label = QGraphicsTextItem(leaf.name or "")
            self.scene.addItem(label)
            label.setPos(max_x + 10, y_leaf - label.boundingRect().height() / 2)

        self.scene.setSceneRect(0, -10, max_x + 150, len(tree.get_terminals()) * 30 + 20)

