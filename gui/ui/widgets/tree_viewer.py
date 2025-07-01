from __future__ import annotations

"""A simple QGraphicsView-based viewer for phylogenetic trees."""

from typing import Callable, Dict, Optional, List, Tuple
import os
from PyQt6.QtWidgets import (
    QWidget,
    QGraphicsView,
    QGraphicsScene,
    QVBoxLayout,
    QGraphicsTextItem,
    QLabel,
    QPushButton,
    QFileDialog,
    QMessageBox,
    QMenu,
)
from PyQt6.QtGui import QPainter, QPen, QColor
from PyQt6.QtCore import Qt
from Bio.Phylo.Newick import Tree, Clade


class _ZoomableGraphicsView(QGraphicsView):
    """Graphics view that supports wheel-based zooming."""

    def wheelEvent(self, event):
        factor = 1.25 if event.angleDelta().y() > 0 else 0.8
        self.scale(factor, factor)

    # ------------------------------------------------------------------
    def contextMenuEvent(self, event):
        parent = self.parentWidget()
        if not isinstance(parent, TreeViewer):
            super().contextMenuEvent(event)
            return
        item = self.itemAt(event.pos())
        if isinstance(item, QGraphicsTextItem) and hasattr(item, "species_name"):
            parent._show_label_menu(item, self.mapToGlobal(event.pos()))
        elif isinstance(item, QGraphicsTextItem) and hasattr(item, "pair_index"):
            parent._show_pair_menu(item, self.mapToGlobal(event.pos()))
        else:
            super().contextMenuEvent(event)


class TreeViewer(QWidget):
    """Window displaying a Newick phylogenetic tree."""

    def __init__(
        self,
        tree: Tree,
        phenotypes: Optional[Dict[str, int]] = None,
        *,
        on_pheno_changed: Optional[Callable[[str], None]] = None,
        on_groups_saved: Optional[Callable[[str], None]] = None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Phylogenetic Tree")
        layout = QVBoxLayout(self)

        if phenotypes:
            legend = QLabel(
                "<b>Legend:</b> <span style='color: blue'>Convergent</span> | "
                "<span style='color: red'>Control</span>"
            )
            layout.addWidget(legend)

        self._phenotypes = phenotypes or {}
        self._tree = tree
        self._on_pheno_changed = on_pheno_changed
        self._on_groups_saved = on_groups_saved

        pheno_btn = QPushButton("Load Phenotype File")
        pheno_btn.clicked.connect(self._select_phenotypes)
        layout.addWidget(pheno_btn)

        self.save_btn = QPushButton("Save Species Groups")
        self.save_btn.setEnabled(False)
        self.save_btn.setToolTip("Must have at least two pairs")
        self.save_btn.clicked.connect(self._save_groups)
        layout.addWidget(self.save_btn)

        self.view = _ZoomableGraphicsView()
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.view.setTransformationAnchor(
            QGraphicsView.ViewportAnchor.AnchorUnderMouse
        )

        layout.addWidget(self.view)
        self.scene = QGraphicsScene(self)
        self.view.setScene(self.scene)

        # pair tracking
        self._pairs: List[Tuple[str, str]] = []
        self._current_role: str | None = None
        self._current_first: str | None = None
        self._disabled_species: set[str] = set()
        self._label_items: Dict[str, QGraphicsTextItem] = {}
        self._branch_lines: Dict[Tuple[Clade, Clade], List] = {}
        self._pair_labels: List[QGraphicsTextItem] = []

        self._draw_tree(tree)

        # Initial window size
        self.resize(1200, 1200)
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

    # ------------------------------------------------------------------
    def _show_label_menu(self, item: QGraphicsTextItem, pos) -> None:
        name = getattr(item, "species_name", "")
        menu = QMenu()
        if name in self._disabled_species:
            act = menu.addAction("Not a valid option")
            act.setEnabled(False)
            menu.exec(pos)
            return

        idx = len(self._pairs) + (1 if self._current_role is None else 0)
        if self._current_role is None:
            conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
            ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
        elif self._current_role == "convergent":
            conv_act = None
            ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
        else:
            ctrl_act = None
            conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
        action = menu.exec(pos)
        if action is None:
            return
        if action == conv_act:
            self._add_species(name, "convergent")
        elif action == ctrl_act:
            self._add_species(name, "control")

    # ------------------------------------------------------------------
    def _show_pair_menu(self, item: QGraphicsTextItem, pos) -> None:
        idx = getattr(item, "pair_index", -1)
        if idx < 1 or idx > len(self._pairs):
            return
        menu = QMenu()
        remove = menu.addAction("Remove Pair")
        action = menu.exec(pos)
        if action == remove:
            self._remove_pair(idx)

    # ------------------------------------------------------------------
    def _update_save_btn(self) -> None:
        self.save_btn.setEnabled(len(self._pairs) >= 2)

    # ------------------------------------------------------------------
    def _add_species(self, name: str, role: str) -> None:
        if self._current_role is None:
            self._current_role = role
            self._current_first = name
        else:
            if role == self._current_role:
                return
            if self._current_role == "convergent":
                conv, ctrl = self._current_first, name
            else:
                conv, ctrl = name, self._current_first
            self._pairs.append((conv, ctrl))
            self._current_role = None
            self._current_first = None
            self._apply_pairs()

    # ------------------------------------------------------------------
    def _path_to(self, ancestor: Clade, leaf: Clade) -> List[Tuple[Clade, Clade]]:
        path = []
        cur = leaf
        while cur is not ancestor:
            parent = self._parent_map.get(cur)
            if parent is None:
                break
            path.append((parent, cur))
            cur = parent
        return path

    # ------------------------------------------------------------------
    def _apply_pairs(self) -> None:
        # reset visuals
        for lines in self._branch_lines.values():
            for l in lines:
                l.setPen(QPen(Qt.GlobalColor.black))
        for name, label in self._label_items.items():
            pheno = self._phenotypes.get(name)
            if pheno == 1:
                color = QColor("blue")
            elif pheno == -1:
                color = QColor("red")
            else:
                color = QColor("black")
            label.setDefaultTextColor(color)
            label.setToolTip("")
        for p in self._pair_labels:
            self.scene.removeItem(p)
        self._pair_labels.clear()
        self._disabled_species.clear()

        # apply all pairs sequentially
        for idx, (conv_name, ctrl_name) in enumerate(self._pairs, start=1):
            conv_leaf = next(self._tree.find_clades(name=conv_name))
            ctrl_leaf = next(self._tree.find_clades(name=ctrl_name))
            ancestor = self._tree.common_ancestor(conv_leaf, ctrl_leaf)

            conv_path = self._path_to(ancestor, conv_leaf)
            ctrl_path = self._path_to(ancestor, ctrl_leaf)
            for parent, child in conv_path:
                for l in self._branch_lines.get((parent, child), []):
                    l.setPen(QPen(QColor("blue"), 2))
            for parent, child in ctrl_path:
                for l in self._branch_lines.get((parent, child), []):
                    l.setPen(QPen(QColor("red"), 2))

            # gray out other descendants
            for leaf in ancestor.get_terminals():
                lname = leaf.name or ""
                if lname not in (conv_name, ctrl_name):
                    lbl = self._label_items.get(lname)
                    if lbl:
                        lbl.setDefaultTextColor(QColor("gray"))
                        lbl.setToolTip("Not a valid option")
                    self._disabled_species.add(lname)

            # label the pair
            x = self._node_pos.get(ancestor, (0, 0))[0]
            y = self._node_pos.get(ancestor, (0, 0))[1]
            label = QGraphicsTextItem(f"Pair {idx}")
            label.pair_index = idx
            self.scene.addItem(label)
            label.setPos(x - label.boundingRect().width() - 5, y)
            self._pair_labels.append(label)

        self._update_save_btn()

    # ------------------------------------------------------------------
    def _remove_pair(self, idx: int) -> None:
        if idx < 1 or idx > len(self._pairs):
            return
        self._pairs.pop(idx - 1)
        self._apply_pairs()

    # ------------------------------------------------------------------
    def _save_groups(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Species Groups",
            os.getcwd(),
            "Text Files (*.txt);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path, "w") as f:
                for conv, ctrl in self._pairs:
                    f.write(f"{conv}\n{ctrl}\n")
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to save file:\n{exc}",
            )
            return
        if hasattr(self, "_on_groups_saved") and self._on_groups_saved:
            self._on_groups_saved(path)


    # ------------------------------------------------------------------
    def _select_phenotypes(self) -> None:
        """Prompt the user for a phenotype file and redraw the tree."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Phenotype File",
            os.getcwd(),
            "CSV Files (*.csv *.txt);;All Files (*)",
        )
        if not path:
            return
        phenos: Dict[str, int] = {}
        try:
            import csv

            with open(path, newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) >= 2:
                        try:
                            phenos[row[0].strip()] = int(row[1])
                        except ValueError:
                            continue
        except Exception as exc:
            QMessageBox.warning(
                self,
                "Phenotypes Error",
                f"Failed to parse phenotypes file:\n{exc}",
            )
            return

        self._phenotypes = phenos
        self.scene.clear()
        self._draw_tree(self._tree)
        if self._on_pheno_changed:
            self._on_pheno_changed(path)

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
        depths = tree.depths()  # parse branch lengths
        max_x = max(depths.values()) if depths else 0
        y_pos = self._y_positions(tree)

        # scale so the tree fills the window horizontally
        page_width = 1200
        label_margin = 150
        px_per_unit = (page_width - label_margin) / max_x if max_x else 1

        def scaled_x(clade: Clade) -> float:
            return depths.get(clade, 0) * px_per_unit

        pen = QPen(Qt.GlobalColor.black)
        self._parent_map = {}
        self._node_pos = {}
        for clade in tree.find_clades(order="preorder"):
            x_parent, y_parent = scaled_x(clade), y_pos.get(clade, 0)
            self._node_pos[clade] = (x_parent, y_parent)
            for child in clade.clades:
                x_child, y_child = scaled_x(child), y_pos.get(child, 0)
                self._node_pos[child] = (x_child, y_child)
                # horizontal segment from parent x to child x at child's y
                h = self.scene.addLine(x_parent, y_child, x_child, y_child, pen)
                # vertical segment from parent y to child y at parent x
                v = self.scene.addLine(x_parent, y_parent, x_parent, y_child, pen)
                self._branch_lines[(clade, child)] = [h, v]
                self._parent_map[child] = clade

        x_max_scaled = max_x * px_per_unit

        # horizontal lines to a common x for tips
        for leaf in tree.get_terminals():
            x_leaf = scaled_x(leaf)
            y_leaf = y_pos.get(leaf, 0)
            line = self.scene.addLine(x_leaf, y_leaf, x_max_scaled, y_leaf, pen)
            self._branch_lines[(leaf, None)] = [line]
            self._node_pos[leaf] = (x_leaf, y_leaf)
            label = QGraphicsTextItem(leaf.name or "")
            label.species_name = leaf.name or ""
            pheno = self._phenotypes.get(leaf.name)
            if pheno == 1:
                label.setDefaultTextColor(QColor("blue"))
            elif pheno == -1:
                label.setDefaultTextColor(QColor("red"))
            self.scene.addItem(label)
            label.setPos(x_max_scaled + 10, y_leaf - label.boundingRect().height() / 2)
            self._label_items[label.species_name] = label

        # Set scene rect based on all items so panning works when zoomed
        bounds = self.scene.itemsBoundingRect()
        self.scene.setSceneRect(bounds.adjusted(-10, -10, 10, 10))
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

