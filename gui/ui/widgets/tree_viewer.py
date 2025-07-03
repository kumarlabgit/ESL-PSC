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
    QGraphicsLineItem,
    QLabel,
    QPushButton,
    QFileDialog,
    QMessageBox,
    QMenu,
    QSizePolicy,
)
from PyQt6.QtGui import QPainter, QPen, QColor, QPalette, QPixmap, QCursor
from PyQt6.QtSvg import QSvgGenerator
from PyQt6.QtCore import Qt, QEvent
from Bio.Phylo.Newick import Tree, Clade


class _ZoomableGraphicsView(QGraphicsView):
    """Graphics view that supports wheel-based zooming with limits."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._zoom_level = 0
        self._min_zoom = -5
        self._max_zoom = 10

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
        if not isinstance(parent, TreeViewer):
            super().contextMenuEvent(event)
            return
        items = self.items(event.pos())
        label_item = None
        pair_item = None
        for it in items:
            if isinstance(it, QGraphicsTextItem) and hasattr(it, "species_name"):
                label_item = it
                break
            if isinstance(it, QGraphicsTextItem) and hasattr(it, "pair_index"):
                pair_item = it
        if label_item is not None:
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            label_item.setCursor(Qt.CursorShape.OpenHandCursor)
            parent._show_label_menu(label_item, self.mapToGlobal(event.pos()))
        elif pair_item is not None:
            self.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            pair_item.setCursor(Qt.CursorShape.OpenHandCursor)
            parent._show_pair_menu(pair_item, self.mapToGlobal(event.pos()))
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
        parent = None
        views = self.scene().views() if self.scene() else []
        if views:
            parent = views[0].parentWidget()
        if isinstance(parent, TreeViewer) and parent._assign_mode in (1, -1):
            cur = parent._assign_cursors.get(parent._assign_mode)
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
        parent = None
        views = self.scene().views() if self.scene() else []
        if views:
            parent = views[0].parentWidget()
        if (
            event.button() == Qt.MouseButton.LeftButton
            and isinstance(parent, TreeViewer)
            and hasattr(self, "species_name")
        ):
            mode = parent._assign_mode
            if mode in (1, -1):
                parent._assign_pheno(self.species_name, mode)
                event.accept()
                return
        super().mousePressEvent(event)


@dataclass
class PairInfo:
    """Information about a convergent-control pair with alternates."""

    convergent: str
    control: str
    conv_alts: List[str] = field(default_factory=list)
    ctrl_alts: List[str] = field(default_factory=list)


@dataclass
class CandidatePair:
    """Helper structure for auto-selected pairs."""

    convergent: str
    control: str
    distance: float
    descendants: set[str] = field(default_factory=set)


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
        self.setWindowTitle("Contrast Pair Selection")
        layout = QVBoxLayout(self)

        # Create a horizontal layout for the top row
        top_row = QHBoxLayout()
        
        # "Set" label and phenotype toggle buttons shown on the right
        set_lbl = QLabel("<b>Set:</b>")
        self.conv_mode_btn = QPushButton("Convergent")
        self.conv_mode_btn.setStyleSheet("color: blue;")
        self.ctrl_mode_btn = QPushButton("Non-convergent")
        self.ctrl_mode_btn.setStyleSheet("color: red;")
        for btn in (self.conv_mode_btn, self.ctrl_mode_btn):
            btn.setCheckable(True)
            btn.setSizePolicy(
                QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
            )
        self.conv_mode_btn.clicked.connect(
            lambda: self._toggle_assign_mode(1)
        )
        self.ctrl_mode_btn.clicked.connect(
            lambda: self._toggle_assign_mode(-1)
        )

        # We will add phenotype buttons first, then a stretch, then the set controls so they stay on the far right
        # Add the top row to the main layout at the end of the setup below (after buttons).
        self._top_row = top_row  # save if needed elsewhere

        self._phenotypes = phenotypes or {}
        self._tree = tree
        self._on_pheno_changed = on_pheno_changed
        self._on_groups_saved = on_groups_saved

        # Branch lines storage must be defined before we compute pen
        self._branch_lines: Dict[Tuple[Clade, Clade | None], List[QGraphicsLineItem]] = {}

        # Phenotype assignment mode tracking and cursors
        self._assign_mode: int | None = None
        self._assign_cursors = {
            1: self._make_color_cursor(QColor("blue")),
            -1: self._make_color_cursor(QColor("red")),
        }


        # Create phenotype buttons
        pheno_btn = QPushButton("Load Phenotype File")
        pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        pheno_btn.setToolTip(
            "Load a CSV file that maps species names (as they appear in the tree) "
            "to trait values: 1 for convergent and -1 for ancestral"
        )
        pheno_btn.clicked.connect(self._select_phenotypes)

        save_pheno_btn = QPushButton("Save Phenotype File")
        save_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        save_pheno_btn.setToolTip("Write the current phenotype assignments to disk")
        save_pheno_btn.clicked.connect(self._save_phenotypes)

        invert_pheno_btn = QPushButton("Invert Phenotype")
        invert_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        invert_pheno_btn.setToolTip(
            "Swap convergent and control phenotypes for all species"
        )
        invert_pheno_btn.clicked.connect(self._invert_phenotypes)

        # Additional phenotype utility buttons
        set_nonconv_btn = QPushButton("Set All to Non-convergent")
        set_nonconv_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        set_nonconv_btn.setToolTip(
            "Assign the Non-convergent phenotype (-1) to every species in the tree"
        )
        set_nonconv_btn.clicked.connect(self._set_all_non_convergent)

        clear_pheno_btn = QPushButton("Clear Phenotypes")
        clear_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        clear_pheno_btn.setToolTip("Remove all phenotype assignments")
        clear_pheno_btn.clicked.connect(self._clear_phenotypes)

        # Add phenotype buttons to the top row
        top_row.addWidget(pheno_btn)
        top_row.addWidget(save_pheno_btn)
        top_row.addWidget(invert_pheno_btn)
        top_row.addWidget(set_nonconv_btn)
        top_row.addWidget(clear_pheno_btn)
        # Stretch so the set buttons stay at far right
        top_row.addStretch()
        top_row.addWidget(set_lbl)
        top_row.addWidget(self.conv_mode_btn)
        top_row.addWidget(self.ctrl_mode_btn)

        # Finally, add this top row to the main layout
        layout.addLayout(top_row)

        export_btn = QPushButton("Export Tree Image")
        export_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        export_btn.setToolTip("Save the current tree view as an SVG graphic")
        export_btn.clicked.connect(self._export_svg)

        self.auto_btn = QPushButton("Auto Select Contrast Pairs")
        self.auto_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.auto_btn.clicked.connect(self._auto_select_pairs)
        self.auto_btn.setEnabled(False)
        self.auto_btn.setToolTip(
            "Automatically choose a valid set of contrast pairs based on the current phenotype assignments"
        )

        groups_btn = QPushButton("Load Species Groups")
        groups_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        groups_btn.setToolTip("Load a saved set of contrast pairs from disk")
        groups_btn.clicked.connect(self._load_groups)

        self.save_btn = QPushButton("Save Species Groups")
        self.save_btn.setEnabled(False)
        self.save_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.save_btn.setToolTip(
            "Save the currently defined convergent/control pairs to a text file"
        )
        self.save_btn.clicked.connect(self._save_groups)

        # ---------------------------
        # Second row: instruction label on left, export & group buttons on right
        bottom_row = QHBoxLayout()
        bottom_row.setContentsMargins(0, 0, 0, 0)
        bottom_row.setSpacing(5)

        instruct = QLabel("Right click species names to add them to the analysis")
        bottom_row.addWidget(instruct)
        bottom_row.addStretch()
        bottom_row.addWidget(self.auto_btn)
        bottom_row.addWidget(export_btn)
        bottom_row.addWidget(groups_btn)
        bottom_row.addWidget(self.save_btn)
        layout.addLayout(bottom_row)

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

        self._pair_labels: List[QGraphicsTextItem] = []
        self._selection_rect: QGraphicsRectItem | None = None
        self._alt_lines: List[Tuple[List, List]] = []
        self._alt_boxes: List[QGraphicsRectItem] = []
        self._main_boxes: Dict[str, QGraphicsRectItem] = {}

        # Determine branch line color based on current palette
        self._update_line_pen()

        # Flag to control whether we should auto-fit the view (only on first draw)
        self._initial_draw = True
        self._draw_tree(tree)

        # Initial window size
        self.resize(1200, 1200)
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def changeEvent(self, event):
        """Update colors when palette (e.g., system light/dark mode) changes."""
        if event.type() == QEvent.Type.PaletteChange:
            self._update_line_pen()
        super().changeEvent(event)

    # ------------------------------------------------------------------
    def _update_line_pen(self):
        """Recalculate branch line pen color based on current palette."""
        pal = self.palette()
        window_lum = pal.color(QPalette.ColorRole.Window).lightness()
        text_lum = pal.color(QPalette.ColorRole.WindowText).lightness()
        self._line_color = (
            Qt.GlobalColor.white if window_lum < text_lum else Qt.GlobalColor.black
        )
        self._line_pen = QPen(self._line_color)
        # Apply new pen to all existing branch lines.
        for items in self._branch_lines.values():
            for it in items:
                it.setPen(self._line_pen)
        # Update text color for unassigned species and pair labels
        default_color = pal.color(QPalette.ColorRole.WindowText)
        for name, lbl in self._label_items.items():
            if self._phenotypes.get(name) not in (1, -1):
                lbl.setDefaultTextColor(default_color)
        for lbl in self._pair_labels:
            lbl.setDefaultTextColor(default_color)

    # ------------------------------------------------------------------
    def _make_color_cursor(self, color: QColor) -> QCursor:
        """Create a small colored-dot cursor for assignment mode."""
        pix = QPixmap(16, 16)
        pix.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pix)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(color)
        painter.setBrush(color)
        painter.drawEllipse(6, 6, 4, 4)
        painter.end()
        return QCursor(pix, 7, 7)

    # ------------------------------------------------------------------
    def _update_assign_cursor(self) -> None:
        """Apply base cursor styles when (de)activating assignment mode."""
        # We keep the view's cursor default (for panning) and reset all label
        # cursors to the open hand.  The colored dot cursor will be shown only
        # when hovering over a label via :class:`_HoverLabelItem`.
        self.view.viewport().unsetCursor()
        for lbl in self._label_items.values():
            lbl.setCursor(Qt.CursorShape.OpenHandCursor)

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
        self.view.verticalScrollBar().setValue(
            self.view.verticalScrollBar().minimum()
        )
        self._update_assign_cursor()

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
            pair = self._pairs[pair_idx - 1]
            make_conv = make_ctrl = remove_alt = None
            if name in pair.conv_alts:
                make_conv = menu.addAction("Make main convergent")
                remove_alt = menu.addAction("Remove from alternates")
            elif name in pair.ctrl_alts:
                make_ctrl = menu.addAction("Make main control")
                remove_alt = menu.addAction("Remove from alternates")
            if remove_alt is not None:
                menu.addSeparator()
            remove_pair = menu.addAction("Remove Pair")
            action = menu.exec(pos)
            if action is None:
                return
            if action == remove_pair:
                self._remove_pair(pair_idx)
            elif action == remove_alt:
                if name in pair.conv_alts:
                    pair.conv_alts.remove(name)
                elif name in pair.ctrl_alts:
                    pair.ctrl_alts.remove(name)
                self._apply_pairs()
            elif action == make_conv:
                if name in pair.conv_alts:
                    pair.conv_alts.remove(name)
                if pair.convergent not in pair.conv_alts:
                    pair.conv_alts.append(pair.convergent)
                pair.convergent = name
                self._apply_pairs()
            elif action == make_ctrl:
                if name in pair.conv_alts:
                    pair.conv_alts.remove(name)
                if name in pair.ctrl_alts:
                    pair.ctrl_alts.remove(name)
                if pair.control not in pair.ctrl_alts:
                    pair.ctrl_alts.append(pair.control)
                pair.control = name
                self._apply_pairs()
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            item.setCursor(Qt.CursorShape.OpenHandCursor)
            self._update_assign_cursor()
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
                if action is None:
                    return
                if action == conv_act:
                    self._add_alternate(name, tgt_idx, "convergent")
                elif action == ctrl_act:
                    self._add_alternate(name, tgt_idx, "control")
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
            item.setCursor(Qt.CursorShape.OpenHandCursor)
            self._update_assign_cursor()
            return

        idx = len(self._pairs) + 1
        pheno = self._phenotypes.get(name)
        allow_conv = pheno != -1
        allow_ctrl = pheno != 1
        conv_act = ctrl_act = cancel_act = None
        if self._current_role is None:
            if allow_conv:
                conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
            if allow_ctrl:
                ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
        elif self._current_role == "convergent":
            if allow_ctrl:
                ctrl_act = menu.addAction(f"Add as control for Pair {idx}")
            if name == self._current_first:
                cancel_act = menu.addAction("Cancel Selection")
        else:
            if allow_conv:
                conv_act = menu.addAction(f"Add as convergent for Pair {idx}")
            if name == self._current_first:
                cancel_act = menu.addAction("Cancel Selection")
        if conv_act is None and ctrl_act is None and cancel_act is None:
            act = menu.addAction("Not a valid option")
            act.setEnabled(False)

        # phenotype change options
        pheno_conv = pheno_ctrl = pheno_clear = None
        menu.addSeparator()
        if pheno == 1:
            pheno_ctrl = menu.addAction("Change to control phenotype")
            pheno_clear = menu.addAction("Remove phenotype")
        elif pheno == -1:
            pheno_conv = menu.addAction("Change to convergent phenotype")
            pheno_clear = menu.addAction("Remove phenotype")
        else:
            pheno_conv = menu.addAction("Set to convergent phenotype")
            pheno_ctrl = menu.addAction("Set to control phenotype")

        action = menu.exec(pos)
        self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
        if action is None:
            self._update_assign_cursor()
            return
        if action == conv_act:
            self._add_species(name, "convergent")
        elif action == ctrl_act:
            self._add_species(name, "control")
        elif action == cancel_act:
            if self._selection_rect is not None:
                self.scene.removeItem(self._selection_rect)
                self._selection_rect = None
            self._current_role = None
            self._current_first = None
            self._update_save_btn()
        elif action == pheno_conv:
            self._phenotypes[name] = 1
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
        elif action == pheno_ctrl:
            self._phenotypes[name] = -1
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
        elif action == pheno_clear:
            self._phenotypes.pop(name, None)
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
        self._update_assign_cursor()

    # ------------------------------------------------------------------
    def _show_pair_menu(self, item: QGraphicsTextItem, pos) -> None:
        idx = getattr(item, "pair_index", -1)
        if idx < 1 or idx > len(self._pairs):
            return
        menu = QMenu()
        remove = menu.addAction("Remove Pair")
        action = menu.exec(pos)
        self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
        item.setCursor(Qt.CursorShape.OpenHandCursor)
        if action == remove:
            self._remove_pair(idx)
        elif action is None:
            self._update_assign_cursor()
            return
        else:
            self._apply_pairs()
        self._update_assign_cursor()

    # ------------------------------------------------------------------
    def _update_save_btn(self) -> None:
        if len(self._pairs) >= 2 and self._current_role is None:
            self.save_btn.setEnabled(True)
            self.save_btn.setToolTip(
                "Save the currently defined convergent/control pairs to a text file"
            )
        else:
            self.save_btn.setEnabled(False)
            if len(self._pairs) < 2:
                self.save_btn.setToolTip("Must have at least two pairs")
            else:
                self.save_btn.setToolTip("There is an unselected pair")

    # ------------------------------------------------------------------
    def _update_auto_btn(self) -> None:
        count = sum(1 for v in self._phenotypes.values() if v in (1, -1))
        if count >= 4:
            self.auto_btn.setEnabled(True)
            self.auto_btn.setToolTip(
                "Automatically choose contrast pairs based on phenotype changes"
            )
        else:
            self.auto_btn.setEnabled(False)
            self.auto_btn.setToolTip("Add more phenotypes")

    # ------------------------------------------------------------------
    def _toggle_assign_mode(self, mode: int) -> None:
        """Toggle phenotype assignment mode using the colored buttons."""
        if self._assign_mode == mode:
            self._assign_mode = None
            self.conv_mode_btn.setChecked(False)
            self.ctrl_mode_btn.setChecked(False)
        else:
            self._assign_mode = mode
            if mode == 1:
                self.conv_mode_btn.setChecked(True)
                self.ctrl_mode_btn.setChecked(False)
            else:
                self.conv_mode_btn.setChecked(False)
                self.ctrl_mode_btn.setChecked(True)
        self._update_assign_cursor()

    # ------------------------------------------------------------------
    def _assign_pheno(self, name: str, val: int) -> None:
        """Set or toggle phenotype for a species and redraw."""
        if val not in (1, -1):
            return
        if self._phenotypes.get(name) == val:
            self._phenotypes.pop(name, None)
        else:
            self._phenotypes[name] = val
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def _add_species(self, name: str, role: str) -> None:
        if self._current_role is None:
            self._current_role = role
            self._current_first = name
            label = self._label_items.get(name)
            if label:
                color = QColor("blue") if role == "convergent" else QColor("red")
                rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                rect.moveTo(label.pos())
                self._selection_rect = self.scene.addRect(rect, QPen(color, 2))
                self._selection_rect.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            self._update_save_btn()
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
                l.setPen(QPen(self._line_color))
                l.setVisible(True)
        for overlay_lines, bases in getattr(self, "_alt_lines", []):
            for item in overlay_lines:
                self.scene.removeItem(item)
            for b in bases:
                b.setVisible(True)
        self._alt_lines.clear()
        for box in getattr(self, "_alt_boxes", []):
            self.scene.removeItem(box)
        self._alt_boxes.clear()
        for box in getattr(self, "_main_boxes", {}).values():
            self.scene.removeItem(box)
        self._main_boxes.clear()
        for name, label in self._label_items.items():
            pheno = self._phenotypes.get(name)
            if pheno == 1:
                color = QColor("blue")
            elif pheno == -1:
                color = QColor("red")
            else:
                # Show unassigned species in the default text color so they
                # remain visible in dark mode.
                color = self.palette().color(QPalette.ColorRole.WindowText)
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
            for l in self._branch_lines.get((conv_leaf, None), []):
                l.setPen(QPen(QColor("blue"), 2))
            for l in self._branch_lines.get((ctrl_leaf, None), []):
                l.setPen(QPen(QColor("red"), 2))

            for name, color in [
                (conv_name, QColor("blue")),
                (ctrl_name, QColor("red")),
            ]:
                label = self._label_items.get(name)
                if label:
                    rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                    rect.moveTo(label.pos())
                    box = self.scene.addRect(rect, QPen(color, 2))
                    box.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                    self._main_boxes[name] = box

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
                    self._disabled_species.add(lname)

            # draw alternate paths
            cc_edges = set(conv_path + ctrl_path)
            for alt_name in pair.conv_alts:
                alt_leaf = next(self._tree.find_clades(name=alt_name))
                full_path = self._path_to(ancestor, alt_leaf)
                trimmed = []
                for edge in full_path:
                    if edge in cc_edges:
                        break
                    trimmed.append(edge)
                for parent, child in trimmed:
                    bases = self._branch_lines.get((parent, child), [])
                    overlay_lines = []
                    for base in bases:
                        base.setVisible(False)
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#87CEFA"), 2, Qt.PenStyle.DashLine),
                        )
                        line.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                        overlay_lines.append(line)
                    self._alt_lines.append((overlay_lines, bases))
                label = self._label_items.get(alt_name)
                if label:
                    rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                    rect.moveTo(label.pos())
                    box = self.scene.addRect(
                        rect, QPen(QColor("#87CEFA"), 2, Qt.PenStyle.DashLine)
                    )
                    box.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                    self._alt_boxes.append(box)

            for alt_name in pair.ctrl_alts:
                alt_leaf = next(self._tree.find_clades(name=alt_name))
                full_path = self._path_to(ancestor, alt_leaf)
                trimmed = []
                for edge in full_path:
                    if edge in cc_edges:
                        break
                    trimmed.append(edge)
                for parent, child in trimmed:
                    bases = self._branch_lines.get((parent, child), [])
                    overlay_lines = []
                    for base in bases:
                        base.setVisible(False)
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#f4aaaa"), 2, Qt.PenStyle.DashLine),
                        )
                        line.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                        overlay_lines.append(line)
                    self._alt_lines.append((overlay_lines, bases))
                label = self._label_items.get(alt_name)
                if label:
                    rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                    rect.moveTo(label.pos())
                    box = self.scene.addRect(
                        rect, QPen(QColor("#f4aaaa"), 2, Qt.PenStyle.DashLine)
                    )
                    box.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                    self._alt_boxes.append(box)

            # label the pair
            x = self._node_pos.get(ancestor, (0, 0))[0]
            y = self._node_pos.get(ancestor, (0, 0))[1]
            label = _HoverLabelItem(f"Pair {idx}")
            label.pair_index = idx
            label.setDefaultTextColor(
                self.palette().color(QPalette.ColorRole.WindowText)
            )
            self.scene.addItem(label)
            label.setPos(x - label.boundingRect().width() - 5, y)
            self._pair_labels.append(label)

        self._update_save_btn()
        self._update_auto_btn()
        self._update_assign_cursor()

    # ------------------------------------------------------------------
    def _remove_pair(self, idx: int) -> None:
        if idx < 1 or idx > len(self._pairs):
            return
        self._pairs.pop(idx - 1)
        self._apply_pairs()

    # ------------------------------------------------------------------
    def _reset_scene(self) -> None:
        """Clear the scene and reset cached graphics items."""
        self.scene.clear()
        self._branch_lines.clear()
        self._label_items.clear()
        self._pair_labels.clear()
        self._alt_lines.clear()
        self._alt_boxes.clear()
        self._main_boxes.clear()
        self._disabled_species.clear()
        self._species_pair_map.clear()
        self._selection_rect = None

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
            if len(lines) < 2 or len(lines) % 2 != 0:
                raise ValueError("File must contain an even number of non-empty lines")
            pairs: List[PairInfo] = []
            for i in range(0, len(lines), 2):
                conv_parts = [s.strip() for s in lines[i].split(',') if s.strip()]
                ctrl_parts = [s.strip() for s in lines[i + 1].split(',') if s.strip()]
                if not conv_parts or not ctrl_parts:
                    raise ValueError(f"Invalid pair entry near lines {i+1}-{i+2}")
                pairs.append(
                    PairInfo(
                        conv_parts[0],
                        ctrl_parts[0],
                        conv_parts[1:],
                        ctrl_parts[1:],
                    )
                )
            if len(pairs) < 2:
                raise ValueError("Need at least two pairs")
        except Exception as exc:
            # Read first 100 characters for context, if possible
            preview = ""
            try:
                with open(path, "r", errors="ignore") as _pf:
                    preview = _pf.read(100)
            except Exception:
                preview = "<unable to read preview>"

            QMessageBox.critical(
                self,
                "Groups Error",
                (
                    "The file format doesn't look right for a species groups file and it couldn't be loaded. "
                    "Are you sure you meant to load this?\n\n"
                    f"File: {os.path.basename(path)}\nFirst 100 characters:\n\n{preview}\n\nDetails: {exc}"
                ),
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
    def _export_svg(self) -> None:
        """Export the current tree view to an SVG file."""
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Tree Image",
            os.getcwd(),
            "SVG Files (*.svg);;All Files (*)",
        )
        if not path:
            return
        try:
            bounds = self.scene.itemsBoundingRect()
            generator = QSvgGenerator()
            generator.setFileName(path)
            generator.setSize(bounds.size().toSize())
            generator.setViewBox(bounds)
            # When exporting in dark mode, temporary switch default-colored
            # elements to black so the SVG looks like light mode output.
            saved_lines = {}
            for lines in self._branch_lines.values():
                for ln in lines:
                    if ln.pen().color() == self._line_color:
                        saved_lines[ln] = ln.pen()
                        ln.setPen(QPen(Qt.GlobalColor.black))

            saved_labels = {}
            for name, lbl in self._label_items.items():
                if self._phenotypes.get(name) not in (1, -1):
                    saved_labels[lbl] = lbl.defaultTextColor()
                    lbl.setDefaultTextColor(QColor("black"))
            for lbl in self._pair_labels:
                saved_labels[lbl] = lbl.defaultTextColor()
                lbl.setDefaultTextColor(QColor("black"))

            painter = QPainter(generator)
            self.scene.render(painter)
            painter.end()

            for item, pen in saved_lines.items():
                item.setPen(pen)
            for lbl, col in saved_labels.items():
                lbl.setDefaultTextColor(col)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to export SVG:\n{exc}",
            )


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
                for idx, row in enumerate(reader, start=1):
                    if not row:
                        continue
                    if len(row) < 2:
                        raise ValueError(f"Line {idx} must have two columns")
                    val = int(row[1])
                    if val not in (1, -1):
                        raise ValueError(f"Invalid value '{row[1]}' on line {idx}")
                    phenos[row[0].strip()] = val
            if not phenos:
                raise ValueError("No valid phenotype entries found")
        except Exception as exc:
            # Read first 100 characters for context, if possible
            preview = ""
            try:
                with open(path, "r", errors="ignore") as _pf:
                    preview = _pf.read(100)
            except Exception:
                preview = "<unable to read preview>"

            QMessageBox.warning(
                self,
                "Phenotypes Error",
                (
                    "The file format doesn't look right for a phenotype file and it couldn't be loaded. "
                    "Are you sure you meant to load this?\n\n"
                    f"File: {os.path.basename(path)}\nFirst 100 characters:\n{preview}\n\nDetails: {exc}"
                ),
            )
            return

        self._phenotypes = phenos
        self._pairs.clear()
        self._current_role = None
        self._current_first = None
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        if self._on_pheno_changed:
            self._on_pheno_changed(path)

    # ------------------------------------------------------------------
    def _save_phenotypes(self) -> None:
        """Save current phenotypes to a CSV file."""
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Phenotype File",
            os.getcwd(),
            "CSV Files (*.csv *.txt);;All Files (*)",
        )
        if not path:
            return
        try:
            import csv

            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                for leaf in self._tree.get_terminals():
                    name = leaf.name or ""
                    val = self._phenotypes.get(name, 0)
                    writer.writerow([name, val])
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to save phenotype file:\n{exc}",
            )
            return
        if self._on_pheno_changed:
            self._on_pheno_changed(path)

    # ------------------------------------------------------------------
    def _invert_phenotypes(self) -> None:
        """Invert phenotype assignments and swap pair roles."""
        for name, val in list(self._phenotypes.items()):
            if val == 1:
                self._phenotypes[name] = -1
            elif val == -1:
                self._phenotypes[name] = 1
        for pair in self._pairs:
            pair.convergent, pair.control = pair.control, pair.convergent
            pair.conv_alts, pair.ctrl_alts = pair.ctrl_alts, pair.conv_alts
        if self._current_role == "convergent":
            self._current_role = "control"
        elif self._current_role == "control":
            self._current_role = "convergent"
        if self._selection_rect is not None and self._current_first is not None:
            color = QColor("blue") if self._current_role == "convergent" else QColor("red")
            self._selection_rect.setPen(QPen(color, 2))
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()

    # ------------------------------------------------------------------
    def _set_all_non_convergent(self) -> None:
        """Assign the non-convergent phenotype (-1) to every species in the tree."""
        for leaf in self._tree.get_terminals():
            self._phenotypes[leaf.name] = -1
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def _clear_phenotypes(self) -> None:
        """Remove all phenotype assignments from the tree."""
        self._phenotypes.clear()
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def _auto_select_pairs(self) -> None:
        """Automatically choose contrast pairs based on phenotype transitions."""
        existing_desc: set[str] = set()
        for pair in self._pairs:
            conv_leaf = next(self._tree.find_clades(name=pair.convergent))
            ctrl_leaf = next(self._tree.find_clades(name=pair.control))
            anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)
            for leaf in anc.get_terminals():
                existing_desc.add(leaf.name or "")

        candidates: List[CandidatePair] = []
        prev_name = None
        prev_pheno = None
        for leaf in self._tree.get_terminals():
            name = leaf.name or ""
            if name in existing_desc:
                continue
            ph = self._phenotypes.get(name)
            if ph not in (1, -1):
                continue
            if prev_name is not None and ph != prev_pheno:
                if prev_pheno == 1:
                    conv, ctrl = prev_name, name
                else:
                    conv, ctrl = name, prev_name
                conv_leaf = next(self._tree.find_clades(name=conv))
                ctrl_leaf = next(self._tree.find_clades(name=ctrl))
                anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)
                dist = self._tree.distance(conv_leaf, ctrl_leaf)
                desc = {l.name or "" for l in anc.get_terminals()}
                candidates.append(CandidatePair(conv, ctrl, dist, desc))
            prev_name = name
            prev_pheno = ph

        candidates.sort(key=lambda c: c.distance)

        added: List[PairInfo] = []
        invalid: set[str] = set(existing_desc)
        while candidates:
            cand = candidates.pop(0)
            if cand.descendants & invalid:
                continue
            added.append(PairInfo(cand.convergent, cand.control))
            invalid.update(cand.descendants)
            candidates = [c for c in candidates if not (c.descendants & invalid)]

        if len(self._pairs) + len(added) < 2:
            QMessageBox.warning(
                self,
                "Auto Select Error",
                "ESL-PSC requires at least 2 valid contrast pairs and if not all of the species are labeled they may need to label more.",
            )
            return

        if added:
            self._pairs.extend(added)
            self._prune_nested_pairs()
            self._current_role = None
            self._current_first = None
            self._apply_pairs()

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

        pen = QPen(self._line_color)
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
            else:
                label.setDefaultTextColor(
                    self.palette().color(QPalette.ColorRole.WindowText)
                )
            self.scene.addItem(label)
            label.setPos(x_max_scaled + 10, y_leaf - label.boundingRect().height() / 2)
            self._label_items[label.species_name] = label

        # Set scene rect based on all items so panning works when zoomed
        bounds = self.scene.itemsBoundingRect()
        self.scene.setSceneRect(bounds.adjusted(-10, -10, 10, 10))

        # Only auto-fit the view on the initial draw. Afterwards preserve the
        # user's zoom level when the tree is redrawn (e.g., after phenotype
        # changes).
        if getattr(self, "_initial_draw", False):
            self.view.resetTransform()
            self.view.fitInView(
                self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio
            )
            # Mark that the initial draw has completed so subsequent redraws do
            # not override the user's zoom.
            self._initial_draw = False

        # Ensure label colors reflect the current palette
        self._update_line_pen()
        self._update_assign_cursor()

