from __future__ import annotations

"""A simple QGraphicsView-based viewer for phylogenetic trees."""

from typing import Callable, Dict, Optional, List, Tuple
from dataclasses import dataclass, field
import os

import random

from PySide6.QtWidgets import (
    QApplication,

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
    QProgressDialog,
    QDialog,
    QInputDialog,
    QStyle,
)
from PySide6.QtGui import (
    QPainter,
    QPen,
    QColor,
    QPalette,
    QPixmap,
    QCursor,
    QLinearGradient,
    QBrush,
    QIcon,
    QKeySequence,
)
from PySide6.QtSvg import QSvgGenerator
from PySide6.QtCore import Qt, QEvent
from Bio.Phylo.Newick import Tree, Clade

from gui.core.fasta_io import read_fasta
from bisect import bisect_left, bisect_right
from gui.ui.widgets.color_utils import viridis_qcolor, viridis_stops_rgb
from gui.ui.widgets.graphics_view import ZoomableGraphicsView
from gui.ui.widgets.label_items import HoverLabelItem
from gui.ui.widgets.dialogs import (
    AutoSelectOptionsDialog as AutoSelectOptionsDialogExt,
    PhenoThresholdDialog,
)


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
        phenotypes: Optional[Dict[str, float]] = None,
        *,
        on_pheno_changed: Optional[Callable[[str], None]] = None,
        on_groups_saved: Optional[Callable[[str], None]] = None,
        on_alignments_changed: Optional[Callable[[str], None]] = None,

        alignments_dir: str = "",
        initial_groups_file: str = "",
        initial_phenotypes_file: str = "",
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
        self.conv_mode_btn.setStyleSheet(
            "QPushButton { color: blue; } QPushButton:disabled { color: gray; }"
        )
        self.ctrl_mode_btn = QPushButton("Non-convergent")
        self.ctrl_mode_btn.setStyleSheet(
            "QPushButton { color: red; } QPushButton:disabled { color: gray; }"
        )
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
        self._on_alignments_changed = on_alignments_changed
        self._alignments_dir = alignments_dir
        # Optionally preload a species groups file on open (if provided by caller)
        self._initial_groups_file = initial_groups_file or ""
        # Optionally preload a phenotypes file on open (if provided by caller)
        self._initial_phenotypes_file = initial_phenotypes_file or ""

        # Continuous phenotype support
        self._continuous_pheno: bool = False
        self._pheno_min: float = 0.0
        self._pheno_max: float = 1.0
        self._update_pheno_mode_and_range()
        # Remember last used continuous thresholds for the session
        self._last_thresh_lower: float | None = None
        self._last_thresh_upper: float | None = None

        # Map explicit species name to total aligned amino acid sequence length.
        self._seq_lengths: Dict[str, int] = {}
        # If the user chose the "Longest Sequence" option during auto-selection, we
        # annotate species labels with their sequence length for display only.
        self._show_seq_lengths: bool = False

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
            "to trait values. Supports binary (-1/1) and continuous floats for gradient coloring."
        )
        pheno_btn.clicked.connect(self._select_phenotypes)

        save_pheno_btn = QPushButton("Save Phenotype File")
        save_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        save_pheno_btn.setToolTip("Write the current phenotype assignments to disk")
        save_pheno_btn.clicked.connect(self._save_phenotypes)

        self.invert_pheno_btn = QPushButton("Invert Phenotype")
        self.invert_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.invert_pheno_btn.setToolTip(
            "Swap convergent and control phenotypes for all species"
        )
        self.invert_pheno_btn.clicked.connect(self._invert_phenotypes)

        # Additional phenotype utility buttons
        self.set_nonconv_btn = QPushButton("Set All to Non-convergent")
        self.set_nonconv_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.set_nonconv_btn.setToolTip(
            "Assign the Non-convergent phenotype (-1) to every species in the tree"
        )
        self.set_nonconv_btn.clicked.connect(self._set_all_non_convergent)

        clear_pheno_btn = QPushButton("Clear Phenotypes")
        clear_pheno_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        clear_pheno_btn.setToolTip("Remove all phenotype assignments")
        clear_pheno_btn.clicked.connect(self._clear_phenotypes)

        # Add phenotype buttons to the top row
        top_row.addWidget(pheno_btn)
        top_row.addWidget(save_pheno_btn)
        top_row.addWidget(self.invert_pheno_btn)
        top_row.addWidget(self.set_nonconv_btn)
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
        # Second row: instruction label on left, UNDO/REDO, auto-select & export/group buttons on right
        bottom_row = QHBoxLayout()
        bottom_row.setContentsMargins(0, 0, 0, 0)
        bottom_row.setSpacing(5)

        instruct = QLabel("Right click species names to add them to the analysis")
        bottom_row.addWidget(instruct)
        bottom_row.addStretch()
        # Undo/Redo buttons (before Auto Select)
        self.undo_btn = QPushButton()
        # Prefer themed curved undo icon; fall back to standard arrow
        undo_icon = QIcon.fromTheme("edit-undo")
        if undo_icon.isNull():
            undo_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_ArrowBack)
        self.undo_btn.setIcon(undo_icon)
        self.undo_btn.setToolTip("Undo")
        self.undo_btn.setEnabled(False)
        self.undo_btn.clicked.connect(self._undo_action)
        self.undo_btn.setShortcut(QKeySequence(QKeySequence.StandardKey.Undo))

        self.redo_btn = QPushButton()
        # Prefer themed curved redo icon; fall back to standard arrow
        redo_icon = QIcon.fromTheme("edit-redo")
        if redo_icon.isNull():
            redo_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_ArrowForward)
        self.redo_btn.setIcon(redo_icon)
        self.redo_btn.setToolTip("Redo")
        self.redo_btn.setEnabled(False)
        self.redo_btn.clicked.connect(self._redo_action)
        self.redo_btn.setShortcut(QKeySequence(QKeySequence.StandardKey.Redo))

        bottom_row.addWidget(self.undo_btn)
        bottom_row.addWidget(self.redo_btn)
        bottom_row.addWidget(self.auto_btn)
        bottom_row.addWidget(export_btn)
        bottom_row.addWidget(groups_btn)
        bottom_row.addWidget(self.save_btn)
        layout.addLayout(bottom_row)

        self.view = ZoomableGraphicsView()
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
        # Ensure phenotype colors and any pre-existing pair overrides are applied
        # on the first render for visual consistency with subsequent updates.
        self._apply_pairs()

        # Snapshot current tree leaf names to detect later tree changes within the
        # same viewer instance. If leaves change (e.g., a new tree loaded upstream),
        # we'll automatically clear any existing pairs to avoid stale state.
        try:
            self._leaf_names_snapshot = {leaf.name or "" for leaf in self._tree.get_terminals()}
        except Exception:
            self._leaf_names_snapshot = set()

        # Initial window size
        self.resize(1200, 1200)
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
        self._update_auto_btn()
        # Ensure controls reflect current phenotype mode on startup
        self._update_pheno_controls_enabled()
        # Initialize undo/redo stacks
        self._undo_stack: List[Dict] = []
        self._redo_stack: List[Dict] = []
        self._update_undo_redo_btns()

        # If an initial phenotypes file path was provided and exists, load it first
        try:
            if self._initial_phenotypes_file and os.path.exists(self._initial_phenotypes_file):
                # Do not push undo stack or notify upstream on initial load
                self._load_phenotypes_from_path(self._initial_phenotypes_file, push_undo=False, notify=False)
        except Exception:
            # Silently ignore preload errors to avoid blocking the viewer from opening
            pass

        # If an initial groups file path was provided and exists, load it now (after phenotypes)
        try:
            if self._initial_groups_file and os.path.exists(self._initial_groups_file):
                self._load_groups_from_path(self._initial_groups_file)
        except Exception:
            # Silently ignore preload errors to avoid blocking the viewer from opening
            pass

    # ------------------------------------------------------------------
    def _update_pheno_controls_enabled(self) -> None:
        """Enable/disable phenotype-setting controls for continuous mode.

        - When continuous phenotypes are active, manual binary assignment is disabled:
          the Convergent/Non-convergent toggle buttons, "Invert Phenotype", and
          "Set All to Non-convergent" buttons are disabled.
        - When binary, these controls are enabled.
        """
        is_cont = getattr(self, "_continuous_pheno", False)
        self.conv_mode_btn.setEnabled(not is_cont)
        self.ctrl_mode_btn.setEnabled(not is_cont)
        if hasattr(self, "invert_pheno_btn"):
            self.invert_pheno_btn.setEnabled(not is_cont)
            if is_cont:
                self.invert_pheno_btn.setToolTip("Disabled for continuous phenotypes")
            else:
                self.invert_pheno_btn.setToolTip("Swap convergent and control phenotypes for all species")
        if hasattr(self, "set_nonconv_btn"):
            self.set_nonconv_btn.setEnabled(not is_cont)
            if is_cont:
                self.set_nonconv_btn.setToolTip("Disabled for continuous phenotypes")
            else:
                self.set_nonconv_btn.setToolTip("Assign the Non-convergent phenotype (-1) to every species in the tree")

    # ------------------------------------------------------------------
    def changeEvent(self, event):
        """Update colors when palette (e.g., system light/dark mode) changes."""
        if event.type() == QEvent.Type.PaletteChange:
            self._update_line_pen()
            # Re-apply pair overlays to restore highlight colors/visibility
            # after the base branch pens have been updated for the new palette.
            # This ensures colored paths and dashed alternates remain visible.
            self._apply_pairs()
            # If the user had an in-progress selection, recreate the selection box
            if self._current_role and self._current_first:
                label = self._label_items.get(self._current_first)
                if label:
                    color = QColor("blue") if self._current_role == "convergent" else QColor("red")
                    rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                    rect.moveTo(label.pos())
                    self._selection_rect = self.scene.addRect(rect, QPen(color, 2))
                    self._selection_rect.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
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
        # Keep pair labels readable with current palette
        default_color = pal.color(QPalette.ColorRole.WindowText)
        for lbl in self._pair_labels:
            lbl.setDefaultTextColor(default_color)
        # Refresh species label colors: unlabeled -> palette text (white in dark mode),
        # labeled -> phenotype-mapped color (viridis/binary), which should remain stable.
        for name, lbl in self._label_items.items():
            lbl.setDefaultTextColor(self._color_for_species(name))

    # ------------------------------------------------------------------
    def _update_pheno_mode_and_range(self) -> None:
        """Detect continuous phenotypes and prepare stats for color mapping.

        For continuous values, we compute:
        - min/max (kept for compatibility elsewhere)
        - a sorted list of phenotype values for percentile-rank mapping
        """
        vals = [float(v) for v in self._phenotypes.values()]
        if not vals:
            self._continuous_pheno = False
            self._pheno_min, self._pheno_max = 0.0, 1.0
            self._pheno_sorted = []
            # Ensure UI controls are updated even when no phenotypes are present
            if hasattr(self, "_update_pheno_controls_enabled"):
                self._update_pheno_controls_enabled()
            return
        self._pheno_min = min(vals)
        self._pheno_max = max(vals)
        # Continuous if any value is not exactly -1 or 1
        self._continuous_pheno = any(v not in (1.0, -1.0) for v in vals)
        # Maintain a sorted list for percentile-based coloring when continuous
        if self._continuous_pheno:
            self._pheno_sorted = sorted(vals)
        else:
            self._pheno_sorted = []
        # Percentile-based coloring does not use min/max; avoid arbitrary padding.
        # Also refresh tool/button enabled state based on phenotype mode.
        if hasattr(self, "_update_pheno_controls_enabled"):
            self._update_pheno_controls_enabled()

    # ------------------------------------------------------------------
    def _percentile_rank(self, val: float) -> float:
        """Compute the percentile rank of a value in the sorted phenotype list."""
        if not self._pheno_sorted:
            return 0.5
        # Center ties: average of left and right insertion points for smoother mapping
        left = bisect_left(self._pheno_sorted, val)
        right = bisect_right(self._pheno_sorted, val)
        idx = 0.5 * (left + right)
        return idx / len(self._pheno_sorted)

    # ------------------------------------------------------------------
    def _map_value_to_color(self, val: float) -> QColor:
        """Map a phenotype value to a color.
        - Continuous: percentile-based Viridis colormap.
        - Binary: 1 -> blue, -1 -> red.
        """
        if self._continuous_pheno:
            # Percentile rank in [0,1]
            p = self._percentile_rank(val)
            return viridis_qcolor(p)
        # binary fallback
        if val == 1 or val == 1.0:
            return QColor("blue")
        if val == -1 or val == -1.0:
            return QColor("red")
        return self.palette().color(QPalette.ColorRole.WindowText)

    # ------------------------------------------------------------------
    def _color_for_species(self, name: str) -> QColor:
        """Base label color for a species based on phenotype and mode."""
        if name not in self._phenotypes:
            # Species without phenotype use the palette text color (white in dark mode)
            return self.palette().color(QPalette.ColorRole.WindowText)
        val = float(self._phenotypes.get(name, 0.0))
        return self._map_value_to_color(val)

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
        # when hovering over a label via :class:`HoverLabelItem`.
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
    def _add_pheno_actions(self, menu: QMenu, name: str):
        """Add phenotype-related actions to a menu for a given species.

        Returns a tuple of actions: (pheno_conv, pheno_ctrl, pheno_clear, edit_val)
        Each element may be None if not applicable.
        """
        pheno_conv = pheno_ctrl = pheno_clear = edit_val = None
        is_empty = len(self._phenotypes) == 0
        is_cont = getattr(self, "_continuous_pheno", False)
        pheno = self._phenotypes.get(name)
        if is_empty:
            # When no phenotypes are assigned at all, allow both paths
            edit_val = menu.addAction("Edit phenotype value…")
            pheno_conv = menu.addAction("Set to convergent phenotype")
            pheno_ctrl = menu.addAction("Set to control phenotype")
        elif is_cont:
            edit_val = menu.addAction("Edit phenotype value…")
        else:
            if pheno == 1:
                pheno_ctrl = menu.addAction("Change to control phenotype")
                pheno_clear = menu.addAction("Remove phenotype")
            elif pheno == -1:
                pheno_conv = menu.addAction("Change to convergent phenotype")
                pheno_clear = menu.addAction("Remove phenotype")
            else:
                pheno_conv = menu.addAction("Set to convergent phenotype")
                pheno_ctrl = menu.addAction("Set to control phenotype")
        return pheno_conv, pheno_ctrl, pheno_clear, edit_val

    # ------------------------------------------------------------------
    def _handle_pheno_action(self, action, name: str, pheno_conv, pheno_ctrl, pheno_clear, edit_val) -> bool:
        """Handle a phenotype action selection.

        Returns True if the action was handled; False otherwise.
        """
        if action == pheno_conv:
            self._push_undo()
            self._phenotypes[name] = 1
            self._update_pheno_mode_and_range()
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
            return True
        elif action == pheno_ctrl:
            self._push_undo()
            self._phenotypes[name] = -1
            self._update_pheno_mode_and_range()
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
            return True
        elif action == pheno_clear:
            self._push_undo()
            self._phenotypes.pop(name, None)
            self._update_pheno_mode_and_range()
            self._reset_scene()
            self._draw_tree(self._tree)
            self._apply_pairs()
            return True
        elif action == edit_val:
            # Continuous-mode edit: prompt for a new numeric value
            cur = float(self._phenotypes.get(name, 0.0)) if self._phenotypes.get(name) is not None else 0.0
            val, ok = QInputDialog.getDouble(self, "Edit phenotype value", f"Set continuous value for {name}:", cur, -1e12, 1e12, 3)
            if ok:
                self._push_undo()
                self._phenotypes[name] = float(val)
                self._update_pheno_mode_and_range()
                self._reset_scene()
                self._draw_tree(self._tree)
                self._apply_pairs()
                self._update_auto_btn()
            return True
        return False

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
            # Always allow phenotype editing actions as well
            menu.addSeparator()
            pheno_conv, pheno_ctrl, pheno_clear, edit_val = self._add_pheno_actions(menu, name)
            action = menu.exec(pos)
            if action is not None:
                if action == remove_pair:
                    self._push_undo()
                    self._remove_pair(pair_idx)
                elif action == remove_alt:
                    self._push_undo()
                    if name in pair.conv_alts:
                        pair.conv_alts.remove(name)
                    elif name in pair.ctrl_alts:
                        pair.ctrl_alts.remove(name)
                    self._apply_pairs()
                elif action == make_conv:
                    self._push_undo()
                    if name in pair.conv_alts:
                        pair.conv_alts.remove(name)
                    if pair.convergent not in pair.conv_alts:
                        pair.conv_alts.append(pair.convergent)
                    pair.convergent = name
                    self._apply_pairs()
                elif action == make_ctrl:
                    self._push_undo()
                    if name in pair.conv_alts:
                        pair.conv_alts.remove(name)
                    if name in pair.ctrl_alts:
                        pair.ctrl_alts.remove(name)
                    if pair.control not in pair.ctrl_alts:
                        pair.ctrl_alts.append(pair.control)
                    pair.control = name
                    self._apply_pairs()
                else:
                    # Handle phenotype actions if selected
                    handled = self._handle_pheno_action(action, name, pheno_conv, pheno_ctrl, pheno_clear, edit_val)
                    if handled:
                        self._update_auto_btn()
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
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
            # Always allow phenotype editing actions as well
            menu.addSeparator()
            pheno_conv, pheno_ctrl, pheno_clear, edit_val = self._add_pheno_actions(menu, name)
            # If neither alt action nor phenotype actions are available, show disabled note
            if conv_act is None and ctrl_act is None and not any([pheno_conv, pheno_ctrl, pheno_clear, edit_val]):
                act = menu.addAction("Not a valid option")
                act.setEnabled(False)
                menu.exec(pos)
            else:
                action = menu.exec(pos)
                if action is not None:
                    if action == conv_act:
                        self._push_undo()
                        self._add_alternate(name, tgt_idx, "convergent")
                    elif action == ctrl_act:
                        self._push_undo()
                        self._add_alternate(name, tgt_idx, "control")
                    else:
                        handled = self._handle_pheno_action(action, name, pheno_conv, pheno_ctrl, pheno_clear, edit_val)
                        if handled:
                            self._update_auto_btn()
            self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
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
        # phenotype change options (always available regardless of pair/disabled)
        menu.addSeparator()
        pheno_conv, pheno_ctrl, pheno_clear, edit_val = self._add_pheno_actions(menu, name)

        action = menu.exec(pos)
        self.view.viewport().setCursor(Qt.CursorShape.OpenHandCursor)
        if action is None:
            self._update_assign_cursor()
            return
        if action == conv_act:
            self._push_undo()
            self._add_species(name, "convergent")
        elif action == ctrl_act:
            self._push_undo()
            self._add_species(name, "control")
        elif action == cancel_act:
            if self._selection_rect is not None:
                self.scene.removeItem(self._selection_rect)
                self._selection_rect = None
            self._current_role = None
            self._current_first = None
            self._update_save_btn()
        else:
            handled = self._handle_pheno_action(action, name, pheno_conv, pheno_ctrl, pheno_clear, edit_val)
            if handled:
                self._update_auto_btn()
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
        if action == remove:
            self._push_undo()
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
        # When continuous phenotypes are present, enable auto-select as long as
        # there are enough species with values to plausibly form pairs. The
        # threshold dialog will refine groups; validation still occurs later.
        if getattr(self, "_continuous_pheno", False):
            count = len(self._phenotypes)
            if count >= 4:
                self.auto_btn.setEnabled(True)
                self.auto_btn.setToolTip(
                    "Automatically choose contrast pairs based on continuous phenotype thresholds"
                )
            else:
                self.auto_btn.setEnabled(False)
                self.auto_btn.setToolTip("Add more phenotypes")
            return

        # Binary phenotype fallback: need at least two convergents and two controls
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
    def _snapshot_state(self) -> Dict:
        """Capture current phenotype and pair state for undo/redo."""
        pairs_copy = [
            PairInfo(p.convergent, p.control, list(p.conv_alts), list(p.ctrl_alts))
            for p in self._pairs
        ]
        return {
            "phenotypes": dict(self._phenotypes),
            "pairs": pairs_copy,
            "current_role": self._current_role,
            "current_first": self._current_first,
        }

    # ------------------------------------------------------------------
    def _restore_state(self, snap: Dict) -> None:
        self._phenotypes = dict(snap.get("phenotypes", {}))
        pairs = snap.get("pairs", [])
        self._pairs = [
            PairInfo(p.convergent, p.control, list(p.conv_alts), list(p.ctrl_alts))
            for p in pairs
        ]
        # Restore in-progress selection state (for initial selection step)
        self._current_role = snap.get("current_role", None)
        self._current_first = snap.get("current_first", None)
        self._update_pheno_mode_and_range()
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        # Recreate selection rectangle if we had an in-progress selection
        if self._current_role and self._current_first:
            label = self._label_items.get(self._current_first)
            if label:
                color = QColor("blue") if self._current_role == "convergent" else QColor("red")
                rect = label.boundingRect().adjusted(-2, 0, 2, 0)
                rect.moveTo(label.pos())
                self._selection_rect = self.scene.addRect(rect, QPen(color, 2))
                self._selection_rect.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        self._update_auto_btn()
        self._update_save_btn()

    # ------------------------------------------------------------------
    def _push_undo(self) -> None:
        self._undo_stack.append(self._snapshot_state())
        # Clear redo stack on new action
        self._redo_stack.clear()
        # Optional cap
        if len(self._undo_stack) > 50:
            self._undo_stack.pop(0)
        self._update_undo_redo_btns()

    # ------------------------------------------------------------------
    def _update_undo_redo_btns(self) -> None:
        self.undo_btn.setEnabled(bool(self._undo_stack))
        self.redo_btn.setEnabled(bool(self._redo_stack))

    # ------------------------------------------------------------------
    def _undo_action(self) -> None:
        if not self._undo_stack:
            return
        cur = self._snapshot_state()
        snap = self._undo_stack.pop()
        self._redo_stack.append(cur)
        self._restore_state(snap)
        self._update_undo_redo_btns()

    # ------------------------------------------------------------------
    def _redo_action(self) -> None:
        if not self._redo_stack:
            return
        cur = self._snapshot_state()
        snap = self._redo_stack.pop()
        self._undo_stack.append(cur)
        self._restore_state(snap)
        self._update_undo_redo_btns()

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
        self._push_undo()
        if self._phenotypes.get(name) == val:
            self._phenotypes.pop(name, None)
        else:
            self._phenotypes[name] = val
        self._update_pheno_mode_and_range()
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
            # Attempt to add the pair
            self._pairs.append(PairInfo(conv, ctrl))
            # Existing logic will remove pairs that are ancestors of others.
            # If that happens to this newly added pair, the addition "fails".
            self._prune_nested_pairs()
            # If the newly added pair was pruned out, inform the user and keep
            # the current selection active so they can choose a different second.
            if not any(p.convergent == conv and p.control == ctrl for p in self._pairs):
                QMessageBox.information(
                    self,
                    "Cannot Add Pair",
                    (
                        "The MRCA of the two selected species already contains one or more existing pairs as descendants, "
                        "so this pair cannot be added. Please remove the conflicting pair(s) first or choose a different combination."
                    ),
                )
                return
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
        def ancestor_for(pair: PairInfo) -> Clade | None:
            conv, ctrl = pair.convergent, pair.control
            try:
                conv_leaf = next(self._tree.find_clades(name=conv))
                ctrl_leaf = next(self._tree.find_clades(name=ctrl))
                return self._tree.common_ancestor(conv_leaf, ctrl_leaf)
            except Exception:
                return None

        ancestors = [ancestor_for(p) for p in self._pairs]
        keep = []
        for i, anc_i in enumerate(ancestors):
            if anc_i is None:
                # Drop invalid pairs referencing species not present on this tree
                continue
            nested = False
            for j, anc_j in enumerate(ancestors):
                if i == j:
                    continue
                if anc_j is None:
                    continue
                if self._is_descendant(anc_i, anc_j):
                    nested = True
                    break
            if not nested:
                keep.append(self._pairs[i])
        if len(keep) != len(self._pairs):
            self._pairs = keep

    # ------------------------------------------------------------------
    def _purge_invalid_pairs(self) -> bool:
        """Remove any pairs and alternates whose species are not in the current tree.

        Returns True if any changes were made.
        """
        leaf_names = {leaf.name or "" for leaf in self._tree.get_terminals()}
        changed = False
        new_pairs: List[PairInfo] = []
        for p in self._pairs:
            if p.convergent not in leaf_names or p.control not in leaf_names:
                changed = True
                continue
            # Filter alternates to only those present
            conv_alts = [n for n in p.conv_alts if n in leaf_names]
            ctrl_alts = [n for n in p.ctrl_alts if n in leaf_names]
            if conv_alts != p.conv_alts or ctrl_alts != p.ctrl_alts:
                changed = True
            new_pairs.append(PairInfo(p.convergent, p.control, conv_alts, ctrl_alts))
        if changed:
            self._pairs = new_pairs
        return changed

    # ------------------------------------------------------------------
    def _clear_pairs_if_tree_changed(self) -> bool:
        """If the current tree leaf-name set differs from the last snapshot,
        clear all pairs to prevent stale state from another tree.

        Returns True if a reset occurred.
        """
        try:
            current = {leaf.name or "" for leaf in self._tree.get_terminals()}
        except Exception:
            current = set()
        snap = getattr(self, "_leaf_names_snapshot", None)
        if snap is None or current != snap:
            self._pairs.clear()
            self._current_role = None
            self._current_first = None
            self._leaf_names_snapshot = current
            return True
        return False

    # ------------------------------------------------------------------
    def _sort_pairs_by_vertical_position(self) -> None:
        """Sort self._pairs by the Y position of each pair's common ancestor.

        This renumbers pairs so labels appear in vertical order down the page.
        """
        def anc_y(pair: PairInfo) -> float:
            try:
                conv_leaf = next(self._tree.find_clades(name=pair.convergent))
                ctrl_leaf = next(self._tree.find_clades(name=pair.control))
                anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)
                return float(self._node_pos.get(anc, (0.0, 0.0))[1])
            except Exception:
                return 0.0
        self._pairs.sort(key=anc_y)

    # ------------------------------------------------------------------
    def _apply_pairs(self) -> None:
        # If the tree's leaf set has changed since last snapshot, clear pairs
        self._clear_pairs_if_tree_changed()
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
            # Baseline label color from phenotype (gradient or binary)
            label.setDefaultTextColor(self._color_for_species(name))
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
        # Ensure we don't keep stale pairs from another tree
        self._purge_invalid_pairs()
        for idx, pair in enumerate(self._pairs, start=1):
            conv_name, ctrl_name = pair.convergent, pair.control
            conv_leaf = next(self._tree.find_clades(name=conv_name), None)
            ctrl_leaf = next(self._tree.find_clades(name=ctrl_name), None)
            if conv_leaf is None or ctrl_leaf is None:
                # Skip invalid pair entries that reference names not present
                # on the current tree instead of crashing.
                continue
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
                    # Ensure pair label color overrides phenotype color
                    label.setDefaultTextColor(color)
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
                        # Use a faded phenotype-aware color so users can see
                        # which disabled species are potential alternates.
                        if getattr(self, "_continuous_pheno", False):
                            # Keep viridis-based color but reduce alpha
                            base = self._color_for_species(lname)
                            faded = QColor(base)
                            faded.setAlpha(160)
                            lbl.setDefaultTextColor(faded)
                        else:
                            ph = self._phenotypes.get(lname)
                            if ph == 1 or ph == 1.0:
                                # light sky blue (matches alt convergent dash color)
                                lbl.setDefaultTextColor(QColor("#87CEFA"))
                            elif ph == -1 or ph == -1.0:
                                # light pink (matches alt control dash color)
                                lbl.setDefaultTextColor(QColor("#f4aaaa"))
                            else:
                                # No phenotype: keep subdued gray
                                lbl.setDefaultTextColor(QColor("gray"))
                    self._disabled_species.add(lname)

            # draw alternate paths
            cc_edges = set(conv_path + ctrl_path)
            for alt_name in pair.conv_alts:
                alt_leaf = next(self._tree.find_clades(name=alt_name), None)
                if alt_leaf is None:
                    continue
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
                # Also overlay the terminal (leaf) segment so dashes reach the tip
                term_bases = self._branch_lines.get((alt_leaf, None), [])
                if term_bases:
                    term_overlays = []
                    for base in term_bases:
                        base.setVisible(False)
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#87CEFA"), 2, Qt.PenStyle.DashLine),
                        )
                        line.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                        term_overlays.append(line)
                    self._alt_lines.append((term_overlays, term_bases))
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
                alt_leaf = next(self._tree.find_clades(name=alt_name), None)
                if alt_leaf is None:
                    continue
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
                # Also overlay the terminal (leaf) segment so dashes reach the tip
                term_bases = self._branch_lines.get((alt_leaf, None), [])
                if term_bases:
                    term_overlays = []
                    for base in term_bases:
                        base.setVisible(False)
                        line = self.scene.addLine(
                            base.line(),
                            QPen(QColor("#f4aaaa"), 2, Qt.PenStyle.DashLine),
                        )
                        line.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
                        term_overlays.append(line)
                    self._alt_lines.append((term_overlays, term_bases))
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
            label = HoverLabelItem(f"Pair {idx}")
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
        if getattr(self, "_show_seq_lengths", False):
            self._update_seq_length_annotations()

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
    def _load_phenotypes_from_path(self, path: str, *, push_undo: bool = True, notify: bool = True) -> bool:
        """Load phenotypes from a CSV path and redraw.

        The file must contain two columns per row: species name and numeric phenotype value.
        Continuous values are supported. Only exact species-name matches present in the
        current tree are applied; others are ignored. Returns True on success, False on failure.
        """
        if not path:
            return False
        phenos: Dict[str, float] = {}
        try:
            import csv
            with open(path, newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row or len(row) < 2:
                        continue
                    name = (row[0] or "").strip()
                    val_str = (row[1] or "").strip()
                    if not name or not val_str:
                        continue
                    try:
                        val = float(val_str)
                    except ValueError:
                        continue
                    phenos[name] = val
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
            return False

        # Keep only exact matches to the current tree's leaf names to avoid stale/mismatched entries
        leaf_names = {leaf.name or "" for leaf in self._tree.get_terminals()}
        phenos = {n: v for n, v in phenos.items() if n in leaf_names}
        if not phenos:
            QMessageBox.warning(
                self,
                "Phenotypes Warning",
                "No phenotype entries matched species present in the current tree.",
            )
            return False

        if push_undo:
            self._push_undo()
        self._phenotypes = phenos
        self._update_pheno_mode_and_range()
        # Loading phenotypes resets existing pairs to avoid inconsistent states
        self._pairs.clear()
        self._current_role = None
        self._current_first = None
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()
        if notify and self._on_pheno_changed:
            self._on_pheno_changed(path)
        return True

    # ------------------------------------------------------------------
    def _load_groups_from_path(self, path: str) -> bool:
        """Load a species groups file from the given path and apply its pairs.

        Returns True on success, False on failure. On success, this also invokes
        the on_groups_saved callback (to propagate the chosen file path upstream),
        prunes nested pairs, and updates the scene overlays.
        """
        if not path:
            return False
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
            return False

        # Validate that ALL species referenced in the groups file exist in the current tree.
        leaf_names = {leaf.name or "" for leaf in self._tree.get_terminals()}
        referenced: set[str] = set()
        for p in pairs:
            referenced.add(p.convergent)
            referenced.add(p.control)
            referenced.update(p.conv_alts)
            referenced.update(p.ctrl_alts)
        missing = [n for n in referenced if n not in leaf_names]
        if missing:
            QMessageBox.warning(
                self,
                "Groups Warning",
                (
                    "The selected species groups file references species that are not present in the current tree, "
                    "so it was not applied.\n\nMissing (first 10 shown):\n- "
                    + "\n- ".join(missing[:10])
                ),
            )
            return False

        self._pairs = pairs
        self._current_role = None
        self._current_first = None
        # Remove any entries not present in this tree, then prune nested ones
        self._purge_invalid_pairs()
        self._prune_nested_pairs()
        self._apply_pairs()
        if self._on_groups_saved:
            self._on_groups_saved(path)
        return True

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
        self._load_groups_from_path(path)

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
                # Only set to black if species has no phenotype; keep binary/continuous colors as-is
                if name not in self._phenotypes:
                    saved_labels[lbl] = lbl.defaultTextColor()
                    lbl.setDefaultTextColor(QColor("black"))
            for lbl in self._pair_labels:
                saved_labels[lbl] = lbl.defaultTextColor()
                lbl.setDefaultTextColor(QColor("black"))

            painter = QPainter(generator)
            self.scene.render(painter)
            # Draw Viridis legend directly onto the SVG using scene coordinates
            # so it is included in the export. Only when continuous phenotypes
            # are active.
            if getattr(self, "_continuous_pheno", False):
                painter.save()
                margin = 8
                bar_w = 120
                bar_h = 10

                # Place legend at top-left of the scene bounds with a small margin
                bar_x = bounds.left() + margin
                # Title on top, then bar below
                fm = painter.fontMetrics()
                title = "Phenotype"
                title_y = bounds.top() + margin + fm.ascent()
                bar_y = title_y + 2

                text_color = QColor("black")
                painter.setPen(QPen(text_color, 1))
                painter.drawText(bar_x, title_y, title)

                # Horizontal Viridis gradient bar (low -> high)
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
        phenos: Dict[str, float] = {}
        try:
            import csv

            with open(path, newline="") as f:
                reader = csv.reader(f)
                for idx, row in enumerate(reader, start=1):
                    if not row:
                        continue
                    if len(row) < 2:
                        # Skip short or malformed rows (e.g., blank lines)
                        continue
                    name = (row[0] or "").strip()
                    val_str = (row[1] or "").strip()
                    if not name or not val_str:
                        continue
                    try:
                        val = float(val_str)
                    except ValueError:
                        # Skip non-numeric rows (e.g., header)
                        continue
                    phenos[name] = val
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

        # State change: push undo snapshot before applying new phenotypes
        self._push_undo()
        self._phenotypes = phenos
        self._update_pheno_mode_and_range()
        self._pairs.clear()
        self._current_role = None
        self._current_first = None
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()
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
        self._apply_pairs()

    # ------------------------------------------------------------------
    def _set_all_non_convergent(self) -> None:
        """Assign the non-convergent phenotype (-1) to every species in the tree."""
        self._push_undo()
        for leaf in self._tree.get_terminals():
            self._phenotypes[leaf.name] = -1
        # Update mode/range first so the redraw uses binary colors immediately
        self._update_pheno_mode_and_range()
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def _clear_phenotypes(self) -> None:
        """Remove all phenotype assignments from the tree."""
        self._push_undo()
        self._phenotypes.clear()
        self._update_pheno_mode_and_range()
        self._reset_scene()
        self._draw_tree(self._tree)
        self._apply_pairs()
        self._update_auto_btn()

    # ------------------------------------------------------------------
    def _auto_select_pairs(self) -> None:
        """Automatically choose contrast pairs based on phenotype transitions.

        For continuous phenotypes, thresholds are used ONLY for the auto-selection
        algorithm (temporary binary mapping) without changing on-screen coloring or
        the continuous legend.
        """
        # If the tree's leaf set has changed since last snapshot, clear pairs now
        self._clear_pairs_if_tree_changed()
        temp_mapping: Dict[str, int] | None = None
        if self._continuous_pheno:
            dlg_thresh = PhenoThresholdDialog(
                list(self._phenotypes.values()), parent=self
            )
            # Pre-populate with last used thresholds for the session (clamped)
            if self._last_thresh_lower is not None and self._last_thresh_upper is not None:
                vmin, vmax = float(self._pheno_min), float(self._pheno_max)
                low = max(vmin, min(vmax, float(self._last_thresh_lower)))
                up = max(vmin, min(vmax, float(self._last_thresh_upper)))
                # Ensure low <= up; if not, reset to dialog defaults
                if low <= up:
                    dlg_thresh.lower_spin.setValue(low)
                    dlg_thresh.upper_spin.setValue(up)
            else:
                # Default both thresholds to the median of phenotype values at session start
                vals = sorted(float(v) for v in self._phenotypes.values())
                if vals:
                    mid = len(vals) // 2
                    if len(vals) % 2 == 1:
                        median = vals[mid]
                    else:
                        median = 0.5 * (vals[mid - 1] + vals[mid])
                    vmin, vmax = float(self._pheno_min), float(self._pheno_max)
                    med = max(vmin, min(vmax, float(median)))
                    # Allow equality by default: both thresholds start at the median
                    dlg_thresh.lower_spin.setValue(med)
                    dlg_thresh.upper_spin.setValue(med)
            if dlg_thresh.exec() != QDialog.DialogCode.Accepted:
                return
            lower = dlg_thresh.lower_threshold
            upper = dlg_thresh.upper_threshold
            # Allow equality; only error if lower > upper
            if lower > upper:
                QMessageBox.warning(
                    self,
                    "Threshold Error",
                    "Lower threshold must not exceed upper threshold",
                )
                return
            # Remember for this session
            self._last_thresh_lower, self._last_thresh_upper = float(lower), float(upper)
            # Build a temporary binary mapping for the algorithm only
            temp_mapping = {
                name: (1 if val > upper else -1)
                for name, val in self._phenotypes.items()
                if (val > upper) or (val < lower)
            }
            if sum(1 for v in temp_mapping.values() if v in (1, -1)) < 4:
                QMessageBox.warning(
                    self,
                    "Threshold Error",
                    "Not enough species outside the thresholds for auto-selection",
                )
                return

        # Ask user how to resolve ambiguous choices
        dlg = AutoSelectOptionsDialogExt(
            bool(self._alignments_dir), bool(self._continuous_pheno), True, parent=self
        )
        dlg.exec()
        if not dlg.choice:
            return
        method = dlg.choice  # "default", "longest", "shortest", "contrast", "composite", or "random"

        # If longest is requested, ensure we have an alignments directory and
        # collect sequence length info for relevant species.
        if method == "longest":
            if not self._alignments_dir and not self._prompt_alignment_dir():
                return
            # For continuous mode, gather lengths for the thresholded set only
            if self._continuous_pheno and temp_mapping is not None:
                backup = self._phenotypes
                try:
                    self._phenotypes = temp_mapping  # temporary scope for length scan
                    self._ensure_sequence_lengths()
                finally:
                    self._phenotypes = backup
            else:
                self._ensure_sequence_lengths()
            # Optionally show sequence-length annotations next to labels
            self._show_seq_lengths = True
            self._update_seq_length_annotations()

        # Avoid overlapping with existing pairs: exclude current descendants
        existing_desc: set[str] = set()
        for pair in self._pairs:
            # Be robust if a previously added pair references a species not in the
            # currently loaded tree (e.g., after loading a new tree or groups file).
            conv_leaf = next(self._tree.find_clades(name=pair.convergent), None)
            ctrl_leaf = next(self._tree.find_clades(name=pair.control), None)
            if conv_leaf is None or ctrl_leaf is None:
                # Skip invalid pairs instead of crashing; the auto-selector will
                # still proceed using only valid state. We avoid mutating
                # self._pairs here to keep this action non-destructive.
                continue
            anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)
            for leaf in anc.get_terminals():
                existing_desc.add(leaf.name or "")

        # Build candidates from adjacent phenotype transitions across tips
        candidates: List[CandidatePair] = []
        prev_name: str | None = None
        prev_pheno = None
        # Choose phenotype source for algorithm (temporary mapping for continuous)
        pheno_for_algo: Dict[str, int] = temp_mapping if temp_mapping is not None else self._phenotypes  # type: ignore[assignment]
        for leaf in self._tree.get_terminals():
            name = leaf.name or ""
            if name in existing_desc:
                continue
            ph = pheno_for_algo.get(name)
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
                # Robust distance: if branch lengths missing, fall back to node count
                try:
                    dist = float(self._tree.distance(conv_leaf, ctrl_leaf))
                except Exception:
                    conv_path = self._path_to(anc, conv_leaf)
                    ctrl_path = self._path_to(anc, ctrl_leaf)
                    dist = float(len(conv_path) + len(ctrl_path))
                desc = {l.name or "" for l in anc.get_terminals()}
                candidates.append(CandidatePair(conv, ctrl, dist, desc))
            prev_name = name
            prev_pheno = ph

        # Candidate ordering: ALWAYS use the default (shortest-distance-first)
        # strategy to maximize the number of non-overlapping pairs. The chosen
        # method only affects which leaf duo is selected within each ancestor.
        candidates.sort(key=lambda c: c.distance)

        # For composite, compute robust trait- and length-based scores across all possible duos
        if method == "composite":
            # Ensure sequence lengths are available. Prompt for an alignments
            # directory if needed, since composite uses lengths for gating.
            if not self._alignments_dir and not self._prompt_alignment_dir():
                return
            # Gather lengths over the thresholded set when in continuous mode to reduce IO
            if self._continuous_pheno and temp_mapping is not None:
                backup = self._phenotypes
                try:
                    self._phenotypes = temp_mapping  # temporary scope for length scan
                    self._ensure_sequence_lengths()
                finally:
                    self._phenotypes = backup
            else:
                self._ensure_sequence_lengths()
            # Optionally annotate labels with lengths like the 'longest' method
            self._show_seq_lengths = True
            self._update_seq_length_annotations()

            # -----------------------------
            # Global precomputations
            # -----------------------------
            # Robust trait scale S_global using MAD and central 80% range
            def _median(vals: List[float]) -> float:
                n = len(vals)
                if n == 0:
                    return 0.0
                s = sorted(vals)
                m = n // 2
                if n % 2 == 1:
                    return s[m]
                return 0.5 * (s[m - 1] + s[m])

            def _percentile(vals: List[float], p: float) -> float:
                # p in [0,1]; linear interpolation between neighbors
                n = len(vals)
                if n == 0:
                    return 0.0
                s = sorted(vals)
                if n == 1:
                    return s[0]
                pos = (n - 1) * max(0.0, min(1.0, p))
                lo = int(pos)
                hi = lo + 1
                if hi >= n:
                    return s[lo]
                frac = pos - lo
                return s[lo] * (1.0 - frac) + s[hi] * frac

            def _mad(vals: List[float], med: float) -> float:
                if not vals:
                    return 0.0
                devs = [abs(x - med) for x in vals]
                return _median(devs)

            trait_vals = []
            try:
                # Use continuous values across all tips with phenotypes
                trait_vals = [float(v) for v in self._phenotypes.values() if v is not None]
            except Exception:
                trait_vals = []

            med_trait = _median(trait_vals)
            mad = _mad(trait_vals, med_trait)
            p10 = _percentile(trait_vals, 0.10)
            p90 = _percentile(trait_vals, 0.90)
            central80 = p90 - p10
            # Fallback floor proportional to the median magnitude
            floor = max(1e-12, 1e-9 * abs(med_trait))
            S_global = max(1.4826 * mad, (central80 / 2.563) if central80 > 0 else 0.0, floor)

            # Median alignment length across all tips with known lengths
            lens_all = list(self._seq_lengths.values())
            L_med = _median(lens_all) if lens_all else None

            # Tie-band and length-gate params
            eps_abs = 0.15
            eps_rel = 0.05
            r_ok = 0.90
            r_bad = 0.50
            epsilon = 0.05

            # Helper to compute distance with branch-length validity check
            def _duo_distance(a_leaf: Clade, b_leaf: Clade) -> float:
                try:
                    anc2 = self._tree.common_ancestor(a_leaf, b_leaf)
                except Exception:
                    return float("inf")
                ap = self._path_to(anc2, a_leaf)
                bp = self._path_to(anc2, b_leaf)
                # Node-count distance is edges in both paths
                edge_count = float(len(ap) + len(bp))
                # Patristic if all child.branch_length are valid numbers
                valid = True
                total = 0.0
                for (_p, ch) in ap + bp:
                    bl = getattr(ch, "branch_length", None)
                    if bl is None:
                        valid = False
                        break
                    try:
                        total += float(bl)
                    except Exception:
                        valid = False
                        break
                return total if valid else edge_count

            # For each candidate ancestor, compute per-duo S and select per tie rules
            best_by_cand: Dict[int, Tuple[float, str, str]] = {}
            for idx, c in enumerate(candidates):
                # Collect eligible leaves under this ancestor
                try:
                    conv_leaf = next(self._tree.find_clades(name=c.convergent))
                    ctrl_leaf = next(self._tree.find_clades(name=c.control))
                    anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)
                except Exception:
                    continue

                convs: List[str] = []
                ctrls: List[str] = []
                for leaf in anc.get_terminals():
                    nm = leaf.name or ""
                    ph = pheno_for_algo.get(nm)
                    if ph == 1:
                        convs.append(nm)
                    elif ph == -1:
                        ctrls.append(nm)
                if not convs or not ctrls:
                    continue

                # Compute per-duo metrics
                per_duo: List[Tuple[str, str, float, float, float, float]] = []  # (a,b,S,dist,L_duo_or_0,T)
                S_max = None
                for a in convs:
                    for b in ctrls:
                        # Trait values must be present
                        try:
                            va = float(self._phenotypes.get(a)) if self._phenotypes.get(a) is not None else None
                            vb = float(self._phenotypes.get(b)) if self._phenotypes.get(b) is not None else None
                        except Exception:
                            va = None
                            vb = None
                        if va is None or vb is None:
                            continue  # skip duo if trait missing
                        diff = va - vb
                        if S_global <= 0.0:
                            continue
                        T = diff / S_global

                        # Length gate using harmonic mean; lenient policy if unknown -> epsilon
                        la = self._seq_lengths.get(a)
                        lb = self._seq_lengths.get(b)
                        if la is None or lb is None or not L_med or L_med <= 0:
                            G_len = epsilon
                            L_duo = 0.0  # for tie-breaker (treat unknown as smallest)
                        else:
                            # Harmonic mean
                            try:
                                L_duo = 2.0 / (1.0 / float(la) + 1.0 / float(lb))
                            except Exception:
                                L_duo = 0.0
                            if L_med and L_med > 0:
                                r = float(L_duo) / float(L_med)
                            else:
                                r = 0.0
                            if r >= r_ok:
                                G_len = 1.0
                            elif r <= r_bad:
                                G_len = epsilon
                            else:
                                G_len = epsilon + (1.0 - epsilon) * ((r - r_bad) / (r_ok - r_bad))

                        S_val = T * G_len

                        # Distance for tie-breaking
                        try:
                            a_leaf = next(self._tree.find_clades(name=a))
                            b_leaf = next(self._tree.find_clades(name=b))
                        except Exception:
                            dist = float("inf")
                        else:
                            dist = _duo_distance(a_leaf, b_leaf)

                        per_duo.append((a, b, S_val, dist, L_duo, T))
                        if S_max is None or S_val > S_max:
                            S_max = S_val

                if not per_duo or S_max is None:
                    continue

                # Tie threshold and selection
                tie_threshold = max(S_max - eps_abs, (1.0 - eps_rel) * S_max)
                tie_set = [t for t in per_duo if t[2] >= tie_threshold]
                # Sort per tie rules: (smallest dist, then largest T, then largest L_duo, then lexicographic (a,b))
                tie_set.sort(key=lambda t: (t[3], -t[5], -t[4], t[0], t[1]))
                a_best, b_best, s_best, _d, _L, _T = tie_set[0]
                best_by_cand[idx] = (s_best, a_best, b_best)

            # Persist chosen duos so _resolve_pair can use them for the composite method
            self._composite_duo: Dict[Tuple[str, str], Tuple[str, str]] = {}
            for idx, (sc, a, b) in best_by_cand.items():
                key = (candidates[idx].convergent, candidates[idx].control)
                self._composite_duo[key] = (a, b)

        added: List[PairInfo] = []
        invalid: set[str] = set(existing_desc)
        while candidates:
            cand = candidates.pop(0)
            if cand.descendants & invalid:
                continue
            pair = self._resolve_pair(cand, method, pheno_for_algo)
            added.append(pair)
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
            # State change: push undo before applying new pairs
            self._push_undo()
            self._pairs.extend(added)
            self._prune_nested_pairs()
            # Renumber pairs to appear in vertical order down the page
            self._sort_pairs_by_vertical_position()
            self._current_role = None
            self._current_first = None
            self._apply_pairs()

    # ------------------------------------------------------------------
    def _prompt_alignment_dir(self) -> bool:
        """Ask the user to choose an alignments directory."""
        msg = QMessageBox(self)
        msg.setWindowTitle("Alignment Directory Required")
        msg.setText(
            "To use the longest sequence method, you must select a directory of alignments for the species with assigned phenotypes on the tree"
        )
        choose_btn = msg.addButton(
            "Choose alignment directory", QMessageBox.ButtonRole.AcceptRole
        )
        msg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        msg.exec()
        if msg.clickedButton() != choose_btn:
            return False

        path = QFileDialog.getExistingDirectory(
            self, "Select Alignment Directory", os.getcwd()
        )
        if not path:
            return False

        self._alignments_dir = path
        if callable(self._on_alignments_changed):
            self._on_alignments_changed(path)
        return True

    # ------------------------------------------------------------------
    def _ensure_sequence_lengths(self) -> None:
        """Gather sequence length totals for all assigned species.

        If the *Longest Sequence* auto-select mode is active, the label
        annotations are refreshed after lengths are gathered.
        """
        if not self._alignments_dir:
            return
        needed = [n for n, v in self._phenotypes.items() if v in (1, -1) and n not in self._seq_lengths]
        if not needed:
            return

        files = [f for f in os.listdir(self._alignments_dir) if f.lower().endswith((".fas", ".fasta", ".fa"))]
        progress = QProgressDialog("Scanning alignments...", "Cancel", 0, len(files), self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)  # show immediately even for quick tasks
        progress.show()
        count = 0
        for fname in files:
            count += 1
            path = os.path.join(self._alignments_dir, fname)
            try:
                records = read_fasta(path)  # returns List[Tuple[str, str]]
            except Exception:
                records = []

            if records:
                for rid, seq in records:
                    sp = rid.split()[0]
                    if sp not in needed:
                        continue
                    # Count non-gap, non-whitespace characters ("-" and "." treated as gaps)
                    aa_count = sum(1 for c in seq if c not in ("-", ".", " ", "\n", "\r"))
                    if aa_count:
                        self._seq_lengths[sp] = self._seq_lengths.get(sp, 0) + aa_count

            # Update progress dialog
            progress.setValue(count)
            progress.setLabelText(f"Scanning alignments... {count}/{len(files)}")
            QApplication.processEvents()  # allow UI to update

            if progress.wasCanceled():
                break
        progress.close()

        # Update on-screen label annotations if requested. This keeps the view in
        # sync even when auto-selection isn't immediately followed by a full
        # redraw.
        if getattr(self, "_show_seq_lengths", False):
            self._update_seq_length_annotations()

    # ------------------------------------------------------------------
    def _update_seq_length_annotations(self) -> None:
        """Refresh label text to include phenotype values (when continuous) and
        optionally sequence-length annotations.

        Notes:
        - We never mutate ``species_name``; only the visible text is updated.
        - If both phenotype value (continuous) and sequence length are present,
          we render as: ``Name (value; len: N)`` with a semicolon separator; the
          length uses commas for readability and switches to scientific notation
          at >= 100,000.
        - For binary phenotypes, we do not append the value.
        """
        for base_name, label in self._label_items.items():
            parts: List[str] = []
            # Append continuous phenotype value if available for this species
            if getattr(self, "_continuous_pheno", False):
                ph_val = self._phenotypes.get(base_name)
                if ph_val is not None:
                    try:
                        fval = float(ph_val)
                        # Show a tidy float (e.g., 0.123, 2.5); no trailing zeros
                        sval = f"{fval:.3f}".rstrip('0').rstrip('.')
                        # Compute viridis-mapped color for this value and color
                        # ONLY the value substring via HTML span. The default
                        # text color (set elsewhere) will still control the name
                        # and any other unstyled text (e.g., length or base).
                        p = self._percentile_rank(fval)
                        vcol = viridis_qcolor(p)
                        vhex = vcol.name()  # #RRGGBB
                        parts.append(f"<span style=\"color:{vhex}\">{sval}</span>")
                    except Exception:
                        pass

            # Append sequence length if enabled and available
            if getattr(self, "_show_seq_lengths", False):
                length = self._seq_lengths.get(base_name)
                if length:
                    # Format with commas; switch to scientific notation at >= 100,000
                    if length >= 100_000:
                        len_fmt = f"{float(length):.3e}"
                    else:
                        len_fmt = f"{length:,}"
                    parts.append(f"len: {len_fmt}")

            # Use '; ' between value and length when both present; otherwise just the one
            if len(parts) > 1:
                extras = f"{parts[0]}; {parts[1]}"
            else:
                extras = parts[0] if parts else ""
            text = base_name if not extras else f"{base_name} (" + extras + ")"
            # Use rich text when we have colored spans; otherwise plain text is fine.
            if any(part.startswith("<span") for part in parts):
                label.setHtml(text)
            else:
                label.setPlainText(text)
            # Keep position stable after text change.
            pos = label.pos()
            label.setPos(pos.x(), pos.y())

        # After text updates, expand the scene rect so zoom/pan includes appended values
        bounds = self.scene.itemsBoundingRect()
        self.scene.setSceneRect(bounds.adjusted(-10, -10, 10, 10))

    # ------------------------------------------------------------------
    def _pick_longest(self, names: List[str]) -> str:
        lengths = {n: self._seq_lengths.get(n, 0) for n in names}
        if not lengths:
            return names[0]
        max_len = max(lengths.values())
        # Deterministic tie-break: first occurrence by input order
        for n in names:
            if lengths.get(n, 0) == max_len:
                return n
        # Fallback (should not happen)
        return names[0]

    # ------------------------------------------------------------------
    def _resolve_pair(self, cand: CandidatePair, method: str, pheno_map: Optional[Dict[str, int]] = None) -> PairInfo:
        conv_leaf = next(self._tree.find_clades(name=cand.convergent))
        ctrl_leaf = next(self._tree.find_clades(name=cand.control))
        anc = self._tree.common_ancestor(conv_leaf, ctrl_leaf)

        # Use provided phenotype mapping if present (e.g., thresholded continuous),
        # otherwise fall back to the viewer's current phenotypes.
        mapping = pheno_map if pheno_map is not None else self._phenotypes

        convs: List[str] = []
        ctrls: List[str] = []
        for leaf in anc.get_terminals():
            nm = leaf.name or ""
            ph = mapping.get(nm)
            if ph == 1:
                convs.append(nm)
            elif ph == -1:
                ctrls.append(nm)

        conv_choice = cand.convergent
        ctrl_choice = cand.control
        if method == "longest":
            if len(convs) > 1:
                conv_choice = self._pick_longest(convs)
            if len(ctrls) > 1:
                ctrl_choice = self._pick_longest(ctrls)
        elif method == "composite":
            # Use the duo chosen during composite scoring if available
            key = (cand.convergent, cand.control)
            duo = getattr(self, "_composite_duo", {}).get(key)
            if duo:
                conv_choice, ctrl_choice = duo
        elif method == "contrast" and getattr(self, "_continuous_pheno", False):
            # For continuous traits: pick extremes among eligible descendants
            if len(convs) > 1:
                conv_choice = max(convs, key=lambda n: float(self._phenotypes.get(n, float("-inf"))))
            if len(ctrls) > 1:
                ctrl_choice = min(ctrls, key=lambda n: float(self._phenotypes.get(n, float("inf"))))
        elif method == "random":
            if len(convs) > 1:
                conv_choice = random.choice(convs)
            if len(ctrls) > 1:
                ctrl_choice = random.choice(ctrls)

        return PairInfo(conv_choice, ctrl_choice)

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
            label = HoverLabelItem(leaf.name or "")
            label.species_name = leaf.name or ""
            # Use gradient/binary color mapping; default black when missing
            label.setDefaultTextColor(self._color_for_species(label.species_name))
            self.scene.addItem(label)
            label.setPos(x_max_scaled + 10, y_leaf - label.boundingRect().height() / 2)
            self._label_items[label.species_name] = label

        # Always refresh label text to include phenotype values (continuous)
        # and, if enabled, sequence-length annotations. Then set scene rect
        # so it includes the full width of the updated labels.
        self._update_seq_length_annotations()
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
