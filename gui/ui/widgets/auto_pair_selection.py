from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, TYPE_CHECKING

import random

from PySide6.QtWidgets import QApplication, QMessageBox, QDialog, QProgressDialog
from PySide6.QtCore import Qt
from Bio.Phylo.Newick import Clade

from gui.ui.widgets.dialogs import (
    AutoSelectOptionsDialog as AutoSelectOptionsDialogExt,
    PhenoThresholdDialog,
)

if TYPE_CHECKING:
    from .tree_viewer import TreeViewer


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
    # MRCA clade for the original boundary transition that produced this candidate.
    #
    # NOTE: This is intentionally stored instead of a precomputed set of descendant
    # leaf names. Materializing descendant sets for every candidate can become
    # quadratic in both time and memory for large trees.
    ancestor: Clade | None = None


def auto_select_pairs(viewer: "TreeViewer") -> None:
    """Automatically choose contrast pairs based on phenotype transitions.

    For continuous phenotypes, thresholds are used ONLY for the auto-selection
    algorithm (temporary binary mapping) without changing on-screen coloring or
    the continuous legend.
    """
    # If the tree's leaf set has changed since last snapshot, clear pairs now
    viewer._clear_pairs_if_tree_changed()
    temp_mapping: Dict[str, int] | None = None
    if viewer._continuous_pheno:
        dlg_thresh = PhenoThresholdDialog(
            list(viewer._phenotypes.values()), parent=viewer
        )
        # Pre-populate with last used thresholds for the session (clamped)
        if viewer._last_thresh_lower is not None and viewer._last_thresh_upper is not None:
            vmin, vmax = float(viewer._pheno_min), float(viewer._pheno_max)
            low = max(vmin, min(vmax, float(viewer._last_thresh_lower)))
            up = max(vmin, min(vmax, float(viewer._last_thresh_upper)))
            # Ensure low <= up; if not, reset to dialog defaults
            if low <= up:
                dlg_thresh.lower_spin.setValue(low)
                dlg_thresh.upper_spin.setValue(up)
        else:
            # Default both thresholds to the median of phenotype values at session start
            vals = sorted(float(v) for v in viewer._phenotypes.values())
            if vals:
                mid = len(vals) // 2
                if len(vals) % 2 == 1:
                    median = vals[mid]
                else:
                    median = 0.5 * (vals[mid - 1] + vals[mid])
                vmin, vmax = float(viewer._pheno_min), float(viewer._pheno_max)
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
                viewer,
                "Threshold Error",
                "Lower threshold must not exceed upper threshold",
            )
            return
        # Remember for this session
        viewer._last_thresh_lower, viewer._last_thresh_upper = float(lower), float(upper)
        # Build a temporary binary mapping for the algorithm only
        temp_mapping = {
            name: (1 if val > upper else -1)
            for name, val in viewer._phenotypes.items()
            if (val > upper) or (val < lower)
        }
        if sum(1 for v in temp_mapping.values() if v in (1, -1)) < 4:
            QMessageBox.warning(
                viewer,
                "Threshold Error",
                "Not enough species outside the thresholds for auto-selection",
            )
            return

    # Ask user how to resolve ambiguous choices
    dlg = AutoSelectOptionsDialogExt(
        bool(viewer._alignments_dir), bool(viewer._continuous_pheno), True, parent=viewer
    )
    dlg.exec()
    if not dlg.choice:
        return
    method = dlg.choice  # "default", "longest", "shortest", "contrast", "composite", or "random"
    # Alternates configuration
    try:
        num_alternates = int(getattr(dlg, "num_alternates", 0))
    except Exception:
        num_alternates = 0
    try:
        max_combos = int(getattr(dlg, "max_combinations", 1))
    except Exception:
        max_combos = 1
    if num_alternates <= 0:
        max_combos = 1  # force single combo when no alternates requested

    # If longest is requested, ensure we have an alignments directory and
    # collect sequence length info for relevant species.
    if method == "longest":
        if not viewer._alignments_dir and not viewer._prompt_alignment_dir():
            return
        # For continuous mode, gather lengths for the thresholded set only
        if viewer._continuous_pheno and temp_mapping is not None:
            backup = viewer._phenotypes
            try:
                viewer._phenotypes = temp_mapping  # temporary scope for length scan
                viewer._ensure_sequence_lengths()
            finally:
                viewer._phenotypes = backup
        else:
            viewer._ensure_sequence_lengths()
        # Optionally show sequence-length annotations next to labels
        viewer._show_seq_lengths = True
        viewer._update_seq_length_annotations()

    progress = QProgressDialog("Auto-selecting contrast pairs...", "Cancel", 0, 0, viewer)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.setMinimumDuration(0)  # show immediately for large trees
    progress.setValue(0)
    progress.show()
    QApplication.processEvents()

    # Build candidates from adjacent phenotype transitions across tips
    candidates: List[CandidatePair] = []
    prev_name: str | None = None
    prev_pheno = None
    prev_leaf: Clade | None = None
    # Choose phenotype source for algorithm (temporary mapping for continuous)
    pheno_for_algo: Dict[str, int] = temp_mapping if temp_mapping is not None else viewer._phenotypes  # type: ignore[assignment]
    terminals = list(viewer._tree.get_terminals())
    progress.setLabelText("Scanning tree for phenotype transitions...")
    progress.setRange(0, len(terminals))
    for i, leaf in enumerate(terminals, start=1):
        name = leaf.name or ""
        ph = pheno_for_algo.get(name)
        if ph not in (1, -1):
            if i % 250 == 0:
                progress.setValue(i)
                QApplication.processEvents()
                if progress.wasCanceled():
                    progress.close()
                    return
            continue
        if prev_name is not None and ph != prev_pheno:
            if prev_pheno == 1:
                conv, ctrl = prev_name, name
                conv_leaf, ctrl_leaf = prev_leaf, leaf
            else:
                conv, ctrl = name, prev_name
                conv_leaf, ctrl_leaf = leaf, prev_leaf
            if conv_leaf is None or ctrl_leaf is None:
                prev_name = name
                prev_pheno = ph
                prev_leaf = leaf
                if i % 250 == 0:
                    progress.setValue(i)
                    QApplication.processEvents()
                    if progress.wasCanceled():
                        progress.close()
                        return
                continue
            anc = viewer._tree.common_ancestor(conv_leaf, ctrl_leaf)
            # Robust distance: if branch lengths missing, fall back to node count
            try:
                dist = float(viewer._tree.distance(conv_leaf, ctrl_leaf))
            except Exception:
                conv_path = viewer._path_to(anc, conv_leaf)
                ctrl_path = viewer._path_to(anc, ctrl_leaf)
                dist = float(len(conv_path) + len(ctrl_path))
            candidates.append(CandidatePair(conv, ctrl, dist, ancestor=anc))
        prev_name = name
        prev_pheno = ph
        prev_leaf = leaf
        if i % 250 == 0:
            progress.setValue(i)
            QApplication.processEvents()
            if progress.wasCanceled():
                progress.close()
                return
    progress.setValue(len(terminals))

    # Candidate ordering: ALWAYS use the default (shortest-distance-first)
    # strategy to maximize the number of non-overlapping pairs. The chosen
    # method only affects which leaf duo is selected within each ancestor.
    candidates.sort(key=lambda c: c.distance)

    # For composite, compute robust trait- and length-based scores across all possible duos
    if method == "composite":
        # Ensure sequence lengths are available. Prompt for an alignments
        # directory if needed, since composite uses lengths for gating.
        if not viewer._alignments_dir and not viewer._prompt_alignment_dir():
            progress.close()
            return
        # Alignment scanning shows its own progress dialog; hide this one to
        # avoid stacking dialogs.
        progress.hide()
        # Gather lengths over the thresholded set when in continuous mode to reduce IO
        if viewer._continuous_pheno and temp_mapping is not None:
            backup = viewer._phenotypes
            try:
                viewer._phenotypes = temp_mapping  # temporary scope for length scan
                viewer._ensure_sequence_lengths()
            finally:
                viewer._phenotypes = backup
        else:
            viewer._ensure_sequence_lengths()
        progress.show()
        # Optionally annotate labels with lengths like the 'longest' method
        viewer._show_seq_lengths = True
        viewer._update_seq_length_annotations()

        progress.setLabelText("Scoring composite candidates...")
        progress.setRange(0, len(candidates))
        progress.setValue(0)

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
            trait_vals = [float(v) for v in viewer._phenotypes.values() if v is not None]
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
        lens_all = list(viewer._seq_lengths.values())
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
                anc2 = viewer._tree.common_ancestor(a_leaf, b_leaf)
            except Exception:
                return float("inf")
            ap = viewer._path_to(anc2, a_leaf)
            bp = viewer._path_to(anc2, b_leaf)
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
        # Store all duo scores so we can rank alternates later
        viewer._composite_per_duo: Dict[Tuple[str, str], List[Tuple[str, str, float, float, float, float]]] = {}
        for idx, c in enumerate(candidates):
            if idx % 5 == 0:
                progress.setValue(idx)
                QApplication.processEvents()
                if progress.wasCanceled():
                    progress.close()
                    return
            # Collect eligible leaves under this ancestor
            try:
                anc = c.ancestor
            except Exception:
                anc = None
            if anc is None:
                try:
                    conv_leaf = viewer._leaf(c.convergent)
                    ctrl_leaf = viewer._leaf(c.control)
                    if conv_leaf is None or ctrl_leaf is None:
                        continue
                    anc = viewer._tree.common_ancestor(conv_leaf, ctrl_leaf)
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
                        va = float(viewer._phenotypes.get(a)) if viewer._phenotypes.get(a) is not None else None
                        vb = float(viewer._phenotypes.get(b)) if viewer._phenotypes.get(b) is not None else None
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
                    la = viewer._seq_lengths.get(a)
                    lb = viewer._seq_lengths.get(b)
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
                        a_leaf = viewer._leaf(a)
                        b_leaf = viewer._leaf(b)
                        if a_leaf is None or b_leaf is None:
                            raise RuntimeError("missing leaf")
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
            key = (candidates[idx].convergent, candidates[idx].control)
            viewer._composite_per_duo[key] = per_duo

        # Persist chosen duos so _resolve_pair can use them for the composite method
        viewer._composite_duo: Dict[Tuple[str, str], Tuple[str, str]] = {}
        for idx, (sc, a, b) in best_by_cand.items():
            key = (candidates[idx].convergent, candidates[idx].control)
            viewer._composite_duo[key] = (a, b)

    added: List[PairInfo] = []
    # Keep scored alternates so we can enforce a global max-combos cap by removing worst first
    alt_scores: List[Tuple[int, str, str, float]] = []  # (pair_index_in_added, role, name, score) role in {"conv","ctrl"}
    # Track clades that must not overlap with newly selected candidates. In a
    # rooted tree, overlap between two clades implies one is an ancestor of the
    # other, so we can check ancestry rather than materializing terminal sets.
    blocked: List[Clade] = []
    for p in viewer._pairs:
        try:
            conv_leaf = viewer._leaf(p.convergent)
            ctrl_leaf = viewer._leaf(p.control)
            if conv_leaf is None or ctrl_leaf is None:
                continue
            blocked.append(viewer._tree.common_ancestor(conv_leaf, ctrl_leaf))
        except Exception:
            continue

    def overlaps_blocked(anc: Clade) -> bool:
        for b in blocked:
            if viewer._is_descendant(b, anc) or viewer._is_descendant(anc, b):
                return True
        return False

    progress.setLabelText("Selecting non-overlapping pairs...")
    progress.setRange(0, len(candidates))
    progress.setValue(0)
    for idx, cand in enumerate(candidates, start=1):
        if idx % 10 == 0:
            progress.setValue(idx)
            QApplication.processEvents()
            if progress.wasCanceled():
                progress.close()
                return
        if cand.ancestor is not None and overlaps_blocked(cand.ancestor):
            continue
        pair = viewer._resolve_pair(cand, method, pheno_for_algo)
        # Optionally compute alternates for this pair
        if num_alternates > 0:
            try:
                conv_leaf = viewer._leaf(pair.convergent)
                ctrl_leaf = viewer._leaf(pair.control)
                if conv_leaf is None or ctrl_leaf is None:
                    raise RuntimeError("missing leaf")
                anc = viewer._tree.common_ancestor(conv_leaf, ctrl_leaf)
            except Exception:
                anc = None
            if anc is not None:
                # Collect eligible descendants by phenotype mapping used for candidate building
                convs: List[str] = []
                ctrls: List[str] = []
                for leaf in anc.get_terminals():
                    nm = leaf.name or ""
                    ph = pheno_for_algo.get(nm)
                    if ph == 1 and nm != pair.convergent:
                        convs.append(nm)
                    elif ph == -1 and nm != pair.control:
                        ctrls.append(nm)

                # Scoring helpers (return higher-is-better score)
                def score_conv_alt(name: str) -> float:
                    # Score candidate alternate convergent when paired with fixed primary control
                    other = pair.control
                    if method == "composite":
                        key = (cand.convergent, cand.control)
                        pd = getattr(viewer, "_composite_per_duo", {}).get(key, [])
                        vals = [S for (a, b, S, _d, _L, _T) in pd if a == name and b == other]
                        return max(vals) if vals else float("-inf")
                    if method == "contrast" and getattr(viewer, "_continuous_pheno", False):
                        try:
                            va = float(viewer._phenotypes.get(name))
                            vb = float(viewer._phenotypes.get(other))
                            return abs(va - vb)
                        except Exception:
                            return float("-inf")
                    if method == "longest":
                        la = float(viewer._seq_lengths.get(name, 0))
                        lb = float(viewer._seq_lengths.get(other, 0))
                        if la <= 0 or lb <= 0:
                            return 0.0
                        try:
                            return 2.0 / (1.0 / la + 1.0 / lb)
                        except Exception:
                            return 0.0
                    if method in ("shortest", "default"):
                        d = viewer._pair_distance_names(name, other)
                        return -d  # smaller distance is better
                    if method == "random":
                        return random.random()
                    # Fallback to distance
                    return -viewer._pair_distance_names(name, other)

                def score_ctrl_alt(name: str) -> float:
                    # Score candidate alternate control when paired with fixed primary convergent
                    other = pair.convergent
                    if method == "composite":
                        key = (cand.convergent, cand.control)
                        pd = getattr(viewer, "_composite_per_duo", {}).get(key, [])
                        vals = [S for (a, b, S, _d, _L, _T) in pd if a == other and b == name]
                        return max(vals) if vals else float("-inf")
                    if method == "contrast" and getattr(viewer, "_continuous_pheno", False):
                        try:
                            va = float(viewer._phenotypes.get(other))
                            vb = float(viewer._phenotypes.get(name))
                            return abs(va - vb)
                        except Exception:
                            return float("-inf")
                    if method == "longest":
                        la = float(viewer._seq_lengths.get(other, 0))
                        lb = float(viewer._seq_lengths.get(name, 0))
                        if la <= 0 or lb <= 0:
                            return 0.0
                        try:
                            return 2.0 / (1.0 / la + 1.0 / lb)
                        except Exception:
                            return 0.0
                    if method in ("shortest", "default"):
                        d = viewer._pair_distance_names(other, name)
                        return -d
                    if method == "random":
                        return random.random()
                    return -viewer._pair_distance_names(other, name)

                # Rank and take top-k alternates
                if convs:
                    convs_scored = [(n, score_conv_alt(n)) for n in convs]
                    convs_scored = [t for t in convs_scored if t[1] != float("-inf")]
                    convs_scored.sort(key=lambda t: t[1], reverse=True)
                    chosen = [n for (n, s) in convs_scored[:num_alternates]]
                    pair.conv_alts.extend(chosen)
                    idx_in_added = len(added)  # 0-based temp index for 'added'
                    for n, s in convs_scored[:num_alternates]:
                        alt_scores.append((idx_in_added, "conv", n, s))
                if ctrls:
                    ctrls_scored = [(n, score_ctrl_alt(n)) for n in ctrls]
                    ctrls_scored = [t for t in ctrls_scored if t[1] != float("-inf")]
                    ctrls_scored.sort(key=lambda t: t[1], reverse=True)
                    chosen = [n for (n, s) in ctrls_scored[:num_alternates]]
                    pair.ctrl_alts.extend(chosen)
                    idx_in_added = len(added)
                    for n, s in ctrls_scored[:num_alternates]:
                        alt_scores.append((idx_in_added, "ctrl", n, s))

        added.append(pair)
        if cand.ancestor is not None:
            blocked.append(cand.ancestor)
    progress.setValue(len(candidates))
    progress.close()

    if len(viewer._pairs) + len(added) < 2:
        QMessageBox.warning(
            viewer,
            "Auto Select Error",
            "ESL-PSC requires at least 2 valid contrast pairs and if not all of the species are labeled they may need to label more.",
        )
        return

    if added:
        # If a max combination cap is set, prune the worst alternates globally until under the cap
        if max_combos is not None and max_combos > 0:
            # Work on a shallow copy of pairs for calculating combo counts during pruning
            tmp_pairs = [PairInfo(p.convergent, p.control, list(p.conv_alts), list(p.ctrl_alts)) for p in added]
            # Sort alternates globally from worst to best (ascending score)
            alt_scores_sorted = sorted(alt_scores, key=lambda x: (x[3], x[1], x[2]))  # by score primarily
            # Remove until combos <= max
            # Include existing pairs in the combo count
            base_existing = [PairInfo(p.convergent, p.control, list(p.conv_alts), list(p.ctrl_alts)) for p in viewer._pairs]
            while viewer._compute_combo_count(base_existing + tmp_pairs, cap=max_combos) > max_combos and alt_scores_sorted:
                idx, role, name, _score = alt_scores_sorted.pop(0)
                if idx < 0 or idx >= len(tmp_pairs):
                    continue
                if role == "conv" and name in tmp_pairs[idx].conv_alts:
                    tmp_pairs[idx].conv_alts.remove(name)
                elif role == "ctrl" and name in tmp_pairs[idx].ctrl_alts:
                    tmp_pairs[idx].ctrl_alts.remove(name)
            # Apply pruned alternates back to added list
            for i, p in enumerate(tmp_pairs):
                added[i].conv_alts = p.conv_alts
                added[i].ctrl_alts = p.ctrl_alts

        # State change: push undo before applying new pairs
        viewer._push_undo()
        viewer._pairs.extend(added)
        viewer._prune_nested_pairs()
        # Renumber pairs to appear in vertical order down the page
        viewer._sort_pairs_by_vertical_position()
        viewer._current_role = None
        viewer._current_first = None
        viewer._apply_pairs()

# ------------------------------------------------------------------
