from __future__ import annotations

"""A simple QGraphicsView-based viewer for phylogenetic trees."""

from typing import Callable, Dict, Optional, List, Tuple
from dataclasses import dataclass, field
import os
from PyQt6.QtWidgets import (
    QWidget,
    QGraphicsView,
    QGraphicsScene,
    QVBoxLayout,
    QHBoxLayout,
    QGraphicsTextItem,
    QGraphicsRectItem,
    QLabel,
    QPushButton,
    QFileDialog,
    QMessageBox,
    QMenu,
    QSizePolicy,
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
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
        elif isinstance(item, QGraphicsTextItem) and hasattr(item, "pair_index"):
            parent._show_pair_menu(item, self.mapToGlobal(event.pos()))
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
        else:
            super().contextMenuEvent(event)


class _HoverLabelItem(QGraphicsTextItem):
    """Text item that changes cursor on hover."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.OpenHandCursor)

    # ------------------------------------------------------------------
    def hoverEnterEvent(self, event):
        self.setCursor(Qt.CursorShape.ArrowCursor)
        super().hoverEnterEvent(event)

    # ------------------------------------------------------------------
    def hoverLeaveEvent(self, event):
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        super().hoverLeaveEvent(event)


@dataclass
class PairInfo:
    """Information about a convergent-control pair with alternates."""

    convergent: str
    control: str
    conv_alts: List[str] = field(default_factory=list)
    ctrl_alts: List[str] = field(default_factory=list)


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
                "<span style='color: red'>Ancestral</span>"
            )
            layout.addWidget(legend, alignment=Qt.AlignmentFlag.AlignRight)

        self._phenotypes = phenotypes or {}
        self._tree = tree
        self._on_pheno_changed = on_pheno_changed
        self._on_groups_saved = on_groups_saved

        pheno_btn = QPushButton("Load Phenotype File")
        pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        pheno_btn.clicked.connect(self._select_phenotypes)

        groups_btn = QPushButton("Load Species Groups")
        groups_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        groups_btn.clicked.connect(self._load_groups)

        self.save_btn = QPushButton("Save Species Groups")
        self.save_btn.setEnabled(False)
        self.save_btn.setToolTip("Must have at least two pairs")
        self.save_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.save_btn.clicked.connect(self._save_groups)

        btn_layout = QHBoxLayout()
        btn_layout.addWidget(pheno_btn)
        btn_layout.addWidget(groups_btn)
        btn_layout.addWidget(self.save_btn)
        layout.addLayout(btn_layout)

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
        self._pairs: List[PairInfo] = []
        self._current_role: str | None = None
        self._current_first: str | None = None
        self._disabled_species: set[str] = set()
        self._species_pair_map: Dict[str, int] = {}
        self._label_items: Dict[str, QGraphicsTextItem] = {}
        self._branch_lines: Dict[Tuple[Clade, Clade], List] = {}
        self._pair_labels: List[QGraphicsTextItem] = []
        self._selection_rect: QGraphicsRectItem | None = None
        self._alt_lines: List = []
        self._alt_boxes: List[QGraphicsRectItem] = []

        self._draw_tree(tree)

        # Initial window size
        self.resize(1200, 1200)
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

    # ------------------------------------------------------------------
    def showEvent(self, event):
        super().showEvent(event)
        bounds = self.scene.itemsBoundingRect()
        view_width = self.view.viewport().width()
        if bounds.width() > 0:
            scale = view_width / bounds.width()
            self.view.resetTransform()
            self.view.scale(scale, scale)
            self.view.centerOn(bounds.center())

    # ------------------------------------------------------------------
    def _show_label_menu(self, item: QGraphicsTextItem, pos) -> None:
        name = getattr(item, "species_name", "")
        menu = QMenu()
        pair_idx = next(
            (
                i + 1
                for i, p in enumerate(self._pairs)
                if name in (p.convergent, p.control, *p.conv_alts, *p.ctrl_alts)
            ),
            None,
        )
        if pair_idx is not None:
            remove = menu.addAction("Remove Pair")
            if menu.exec(pos) == remove:
                self._remove_pair(pair_idx)
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            return
        if name in self._disabled_species:
            tgt_idx = self._species_pair_map.get(name)
            pheno = self._phenotypes.get(name)
            conv_act = ctrl_act = None
            if tgt_idx is not None:
                if pheno != -1:
                    conv_act = menu.addAction(
                        f"Add as alternate convergent for Pair {tgt_idx}"
                    )
                if pheno != 1:
                    ctrl_act = menu.addAction(
                        f"Add as alternate control for Pair {tgt_idx}"
                    )
            if conv_act is None and ctrl_act is None:
                act = menu.addAction("Not a valid option")
                act.setEnabled(False)
                menu.exec(pos)
            else:
                action = menu.exec(pos)
                if action == conv_act:
                    self._add_alternate(name, tgt_idx, "convergent")
                elif action == ctrl_act:
                    self._add_alternate(name, tgt_idx, "control")
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            return

        idx = len(self._pairs) + (1 if self._current_role is None else 0)
        pheno = self._phenotypes.get(name)
        allow_conv = pheno != -1
        allow_ctrl = pheno != 1
        conv_act = ctrl_act = None
        if self._current_role is None:
            if allow_conv:
                conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
            if allow_ctrl:
                ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
        elif self._current_role == "convergent":
            if allow_ctrl:
                ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
        else:
            if allow_conv:
                conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
        if conv_act is None and ctrl_act is None:
            act = menu.addAction("Not a valid option")
            act.setEnabled(False)
        action = menu.exec(pos)
        self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
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
        self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
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
            label = self._label_items.get(name)
            if label:
                color = QColor("blue") if role == "convergent" else QColor("red")
                rect = label.boundingRect().adjusted(-2, -2, 2, 2)
                rect.moveTo(label.pos())
                self._selection_rect = self.scene.addRect(rect, QPen(color, 2))
        else:
            if role == self._current_role:
                return
            if self._current_role == "convergent":
                conv, ctrl = self._current_first, name
            else:
                conv, ctrl = name, self._current_first
            self._pairs.append(PairInfo(conv, ctrl))
            self._prune_nested_pairs()
            self._current_role = None
            self._current_first = None
            if self._selection_rect is not None:
                self.scene.removeItem(self._selection_rect)
                self._selection_rect = None
            self._apply_pairs()

    # ------------------------------------------------------------------
    def _add_alternate(self, name: str, pair_idx: int, role: str) -> None:
        if pair_idx < 1 or pair_idx > len(self._pairs):
            return
        pair = self._pairs[pair_idx - 1]
        if role == "convergent":
            if name not in pair.conv_alts:
                pair.conv_alts.append(name)
        else:
            if name not in pair.ctrl_alts:
                pair.ctrl_alts.append(name)
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
    def _is_descendant(self, ancestor: Clade, node: Clade) -> bool:
        while node is not None and node is not ancestor:
            node = self._parent_map.get(node)
        return node is ancestor

    # ------------------------------------------------------------------
    def _prune_nested_pairs(self) -> None:
        """Remove pairs that are nested within another pair's ancestor."""
        def ancestor_for(pair: PairInfo) -> Clade:
            conv, ctrl = pair.convergent, pair.control
            conv_leaf = next(self._tree.find_clades(name=conv))
            ctrl_leaf = next(self._tree.find_clades(name=ctrl))
            return self._tree.common_ancestor(conv_leaf, ctrl_leaf)

        ancestors = [ancestor_for(p) for p in self._pairs]
        keep = []
        for i, anc_i in enumerate(ancestors):
            nested = False
            for j, anc_j in enumerate(ancestors):
                if i == j:
                    continue
                if self._is_descendant(anc_i, anc_j):
                    nested = True
                    break
            if not nested:
                keep.append(self._pairs[i])
        if len(keep) != len(self._pairs):
            self._pairs = keep

    # ------------------------------------------------------------------
    def _apply_pairs(self) -> None:
        # reset visuals
        for lines in self._branch_lines.values():
            for l in lines:
                l.setPen(QPen(Qt.GlobalColor.black))
        for line in getattr(self, "_alt_lines", []):
            self.scene.removeItem(line)
        self._alt_lines.clear()
        for box in getattr(self, "_alt_boxes", []):
            self.scene.removeItem(box)
        self._alt_boxes.clear()
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
        self._species_pair_map.clear()
        if self._selection_rect is not None:
            self.scene.removeItem(self._selection_rect)
            self._selection_rect = None

        # apply all pairs sequentially
        for idx, pair in enumerate(self._pairs, start=1):
            conv_name, ctrl_name = pair.convergent, pair.control
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

            # gray out other descendants but allow alternates
            excluded = {conv_name, ctrl_name, *pair.conv_alts, *pair.ctrl_alts}
            for leaf in ancestor.get_terminals():
                lname = leaf.name or ""
                self._species_pair_map[lname] = idx
                if lname not in excluded:
                    lbl = self._label_items.get(lname)
                    if lbl:
                        ph = self._phenotypes.get(lname)
                        if ph == 1:
                            c = QColor("#add8e6")
                        elif ph == -1:
                            c = QColor("#f4aaaa")
                        else:
                            c = QColor("gray")
                        lbl.setDefaultTextColor(c)
                        lbl.setToolTip("Not a valid option")
                    self._disabled_species.add(lname)

            # draw alternate paths
            for alt_name in pair.conv_alts:
                alt_leaf = next(self._tree.find_clades(name=alt_name))
                alt_path = self._path_to(ancestor, alt_leaf)
                for parent, child in alt_path:
                    for base in self._branch_lines.get((parent, child), []):
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#87CEFA"), 2, Qt.PenStyle.DashLine),
                        )
                        self._alt_lines.append(line)
                label = self._label_items.get(alt_name)
                if label:
                    rect = label.boundingRect().adjusted(-2, -2, 2, 2)
                    rect.moveTo(label.pos())
                    box = self.scene.addRect(
                        rect, QPen(QColor("#87CEFA"), 2, Qt.PenStyle.DashLine)
                    )
                    self._alt_boxes.append(box)

            for alt_name in pair.ctrl_alts:
                alt_leaf = next(self._tree.find_clades(name=alt_name))
                alt_path = self._path_to(ancestor, alt_leaf)
                for parent, child in alt_path:
                    for base in self._branch_lines.get((parent, child), []):
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#f4aaaa"), 2, Qt.PenStyle.DashLine),
                        )
                        self._alt_lines.append(line)
                label = self._label_items.get(alt_name)
                if label:
                    rect = label.boundingRect().adjusted(-2, -2, 2, 2)
                    rect.moveTo(label.pos())
                    box = self.scene.addRect(
                        rect, QPen(QColor("#f4aaaa"), 2, Qt.PenStyle.DashLine)
                    )
                    self._alt_boxes.append(box)

            # label the pair
            x = self._node_pos.get(ancestor, (0, 0))[0]
            y = self._node_pos.get(ancestor, (0, 0))[1]
            label = _HoverLabelItem(f"Pair {idx}")
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
    def _load_groups(self) -> None:
        """Load a species groups file and apply its pairs."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Species Groups",
            os.getcwd(),
            "Text Files (*.txt);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            if len(lines) % 2 != 0:
                raise ValueError("File must have an even number of lines")
            pairs: List[PairInfo] = []
            for i in range(0, len(lines), 2):
                conv_parts = [s.strip() for s in lines[i].split(',') if s.strip()]
                ctrl_parts = [s.strip() for s in lines[i + 1].split(',') if s.strip()]
                if not conv_parts or not ctrl_parts:
                    raise ValueError("Invalid pair entry")
                pairs.append(
                    PairInfo(
                        conv_parts[0],
                        ctrl_parts[0],
                        conv_parts[1:],
                        ctrl_parts[1:],
                    )
                )
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Groups Error",
                f"Failed to load species groups:\n{exc}",
            )
            return

        self._pairs = pairs
        self._current_role = None
        self._current_first = None
        self._prune_nested_pairs()
        self._apply_pairs()
        if self._on_groups_saved:
            self._on_groups_saved(path)

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
                for pair in self._pairs:
                    conv_line = ",".join([pair.convergent] + pair.conv_alts)
                    ctrl_line = ",".join([pair.control] + pair.ctrl_alts)
                    f.write(f"{conv_line}\n{ctrl_line}\n")
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
            label = _HoverLabelItem(leaf.name or "")
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

