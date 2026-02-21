from __future__ import annotations

import csv
import os
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from gui.core.fasta_io import read_fasta


@dataclass
class PairInfo:
    convergent: str
    control: str
    conv_alts: List[str] = field(default_factory=list)
    ctrl_alts: List[str] = field(default_factory=list)


@dataclass
class CandidatePair:
    convergent: str
    control: str
    distance: float
    ancestor: Optional[Clade] = None
    pct_diff: float = 0.0


def load_tree(tree_path: str) -> Tree:
    if not os.path.exists(tree_path):
        raise FileNotFoundError(tree_path)
    ext = os.path.splitext(tree_path)[1].lower()
    if ext in {".nexus", ".nex"}:
        return Phylo.read(tree_path, "nexus")
    with open(tree_path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline().strip()
    if first.upper().startswith("#NEXUS"):
        return Phylo.read(tree_path, "nexus")
    return Phylo.read(tree_path, "newick")


def load_phenotypes(pheno_path: str) -> Dict[str, float]:
    if not os.path.exists(pheno_path):
        raise FileNotFoundError(pheno_path)
    phenos: Dict[str, float] = {}
    header_skipped = False
    with open(pheno_path, newline="", encoding="utf-8", errors="ignore") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if not row:
                continue
            if len(row) < 2:
                continue
            sp = (row[0] or "").strip()
            val_s = (row[1] or "").strip()
            if not sp or not val_s:
                continue
            sp_norm = sp.strip().strip('"').strip("'")
            if sp_norm.lower() == "species":
                header_skipped = True
                continue
            try:
                phenos[sp_norm] = float(val_s)
            except ValueError:
                if not header_skipped:
                    header_skipped = True
                    continue
    if not phenos:
        raise ValueError("No valid phenotype entries found")
    return phenos


def _pheno_is_continuous(values: Sequence[float]) -> bool:
    return any(float(v) not in (-1.0, 0.0, 1.0) for v in values)


def _median(vals: Sequence[float]) -> float:
    s = sorted(float(x) for x in vals)
    n = len(s)
    if n == 0:
        return 0.0
    m = n // 2
    if n % 2 == 1:
        return s[m]
    return 0.5 * (s[m - 1] + s[m])


def _percentile(vals: Sequence[float], p: float) -> float:
    s = sorted(float(x) for x in vals)
    n = len(s)
    if n == 0:
        return 0.0
    if n == 1:
        return s[0]
    pos = (n - 1) * max(0.0, min(1.0, p))
    lo = int(pos)
    hi = lo + 1
    if hi >= n:
        return s[lo]
    frac = pos - lo
    return s[lo] * (1.0 - frac) + s[hi] * frac


def thresholds_from_quantile_tails(values: Sequence[float], tails_pct: float) -> Tuple[float, float]:
    q = max(0.0, min(50.0, float(tails_pct))) / 100.0
    if q <= 0.0:
        return min(values), max(values)
    return _percentile(values, q), _percentile(values, 1.0 - q)


def _build_parent_map(root: Clade) -> Dict[Clade, Clade]:
    parent: Dict[Clade, Clade] = {}

    def rec(node: Clade) -> None:
        for ch in getattr(node, "clades", []) or []:
            parent[ch] = node
            rec(ch)

    rec(root)
    return parent


def _path_to(parent_map: Dict[Clade, Clade], ancestor: Clade, leaf: Clade) -> List[Tuple[Clade, Clade]]:
    path: List[Tuple[Clade, Clade]] = []
    cur = leaf
    while cur is not ancestor:
        p = parent_map.get(cur)
        if p is None:
            break
        path.append((p, cur))
        cur = p
    return path


def _is_descendant(parent_map: Dict[Clade, Clade], ancestor: Clade, node: Clade) -> bool:
    cur = node
    while cur is not None and cur is not ancestor:
        cur = parent_map.get(cur)
    return cur is ancestor


def _pair_distance(tree: Tree, parent_map: Dict[Clade, Clade], a: Clade, b: Clade) -> float:
    try:
        return float(tree.distance(a, b))
    except Exception:
        pass
    try:
        anc = tree.common_ancestor(a, b)
    except Exception:
        return float("inf")
    ap = _path_to(parent_map, anc, a)
    bp = _path_to(parent_map, anc, b)
    return float(len(ap) + len(bp))


def _duo_distance(tree: Tree, parent_map: Dict[Clade, Clade], a: Clade, b: Clade) -> float:
    try:
        anc2 = tree.common_ancestor(a, b)
    except Exception:
        return float("inf")
    ap = _path_to(parent_map, anc2, a)
    bp = _path_to(parent_map, anc2, b)
    edge_count = float(len(ap) + len(bp))
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


def _pair_distance_names(tree: Tree, parent_map: Dict[Clade, Clade], leaf_by_name: Dict[str, Clade], a: str, b: str) -> float:
    la = leaf_by_name.get(a)
    lb = leaf_by_name.get(b)
    if la is None or lb is None:
        return float("inf")
    return _pair_distance(tree, parent_map, la, lb)


def _compute_combo_count(pairs: List[PairInfo], cap: Optional[int] = None) -> int:
    total = 1
    for p in pairs:
        choices = (1 + len(p.conv_alts)) * (1 + len(p.ctrl_alts))
        total *= int(choices)
        if cap is not None and total > cap:
            return total
    return total


def _build_pct_contrast_candidates(
    tree: Tree,
    parent_map: Dict[Clade, Clade],
    leaf_by_name: Dict[str, Clade],
    phenotypes: Dict[str, float],
) -> List[CandidatePair]:
    names = [nm for nm in leaf_by_name.keys() if nm in phenotypes]
    out: List[CandidatePair] = []
    for i in range(len(names)):
        a = names[i]
        va = float(phenotypes[a])
        if va <= 0.0:
            continue
        la = leaf_by_name.get(a)
        if la is None:
            continue
        for j in range(i + 1, len(names)):
            b = names[j]
            vb = float(phenotypes[b])
            if vb <= 0.0:
                continue
            lb = leaf_by_name.get(b)
            if lb is None:
                continue
            if va == vb:
                continue

            if va > vb:
                conv, ctrl = a, b
                upper, lower = va, vb
                conv_leaf, ctrl_leaf = la, lb
            else:
                conv, ctrl = b, a
                upper, lower = vb, va
                conv_leaf, ctrl_leaf = lb, la

            if lower <= 0.0:
                continue
            pct_diff = ((upper - lower) / lower) * 100.0
            try:
                anc = tree.common_ancestor(conv_leaf, ctrl_leaf)
            except Exception:
                anc = None
            dist = _pair_distance(tree, parent_map, conv_leaf, ctrl_leaf)
            out.append(
                CandidatePair(
                    convergent=conv,
                    control=ctrl,
                    distance=float(dist),
                    ancestor=anc,
                    pct_diff=float(pct_diff),
                )
            )
    out.sort(key=lambda c: (c.distance, -c.pct_diff, c.convergent, c.control))
    return out


def _select_pct_contrast_candidates(
    candidates: List[CandidatePair],
    parent_map: Dict[Clade, Clade],
    min_pct_diff: float,
    blocked: Optional[List[Clade]] = None,
) -> List[CandidatePair]:
    selected: List[CandidatePair] = []
    blocked_local: List[Clade] = list(blocked or [])

    def overlaps_blocked(anc: Clade) -> bool:
        for b in blocked_local:
            if _is_descendant(parent_map, b, anc) or _is_descendant(parent_map, anc, b):
                return True
        return False

    for cand in candidates:
        if cand.pct_diff < min_pct_diff:
            continue
        anc = cand.ancestor
        if anc is None:
            continue
        if overlaps_blocked(anc):
            continue
        selected.append(cand)
        blocked_local.append(anc)
    return selected


def sweep_pct_contrast_pair_counts(
    tree: Tree,
    phenotypes: Dict[str, float],
    thresholds: Sequence[float],
) -> List[Tuple[float, int]]:
    leaf_names = {leaf.name or "" for leaf in tree.get_terminals()}
    phenos = {n: float(v) for n, v in phenotypes.items() if n in leaf_names}
    if not phenos:
        raise ValueError("No phenotype entries matched species present in the tree")
    vals = list(phenos.values())
    if not _pheno_is_continuous(vals):
        raise ValueError("Percent-contrast thresholding is only available for continuous phenotypes")
    if any(float(v) <= 0.0 for v in vals):
        raise ValueError("Percent-contrast thresholding requires strictly positive phenotype values")

    parent_map = _build_parent_map(tree.root)
    leaf_by_name: Dict[str, Clade] = {}
    for leaf in tree.get_terminals():
        nm = leaf.name or ""
        if nm and nm not in leaf_by_name:
            leaf_by_name[nm] = leaf

    candidates = _build_pct_contrast_candidates(tree, parent_map, leaf_by_name, phenos)
    out: List[Tuple[float, int]] = []
    for th in thresholds:
        thr = max(0.0, float(th))
        cnt = len(_select_pct_contrast_candidates(candidates, parent_map, thr))
        out.append((thr, int(cnt)))
    return out


def _y_positions(tree: Tree, step: int = 30) -> Dict[Clade, float]:
    y: Dict[Clade, float] = {}
    for idx, leaf in enumerate(tree.get_terminals()):
        y[leaf] = idx * step

    def set_internal(clade: Clade) -> float:
        if clade.is_terminal():
            return y[clade]
        vals = [set_internal(c) for c in clade.clades]
        y[clade] = sum(vals) / len(vals) if vals else 0.0
        return y[clade]

    set_internal(tree.root)
    return y


def _scan_sequence_lengths(alignments_dir: str, needed: List[str]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not alignments_dir or not os.path.isdir(alignments_dir) or not needed:
        return out
    needed_set = set(needed)
    files = [
        f
        for f in os.listdir(alignments_dir)
        if f.lower().endswith((".fas", ".fasta", ".fa", ".faa"))
    ]
    for fname in files:
        path = os.path.join(alignments_dir, fname)
        try:
            records = read_fasta(path)
        except Exception:
            records = []
        for rid, seq in records:
            sp = (rid or "").split()[0]
            if sp not in needed_set:
                continue
            aa_count = sum(1 for c in seq if c not in ("-", ".", " ", "\n", "\r"))
            if aa_count:
                out[sp] = out.get(sp, 0) + aa_count
    return out


def _pick_longest(seq_lengths: Dict[str, int], names: List[str]) -> str:
    lengths = {n: seq_lengths.get(n, 0) for n in names}
    if not lengths:
        return names[0]
    max_len = max(lengths.values())
    for n in names:
        if lengths.get(n, 0) == max_len:
            return n
    return names[0]


def _resolve_pair(
    tree: Tree,
    leaf_by_name: Dict[str, Clade],
    seq_lengths: Dict[str, int],
    cand: CandidatePair,
    method: str,
    pheno_for_algo: Dict[str, int],
    phenotypes: Dict[str, float],
    composite_duo: Dict[Tuple[str, str], Tuple[str, str]],
    continuous: bool,
) -> PairInfo:
    conv_leaf = leaf_by_name.get(cand.convergent)
    ctrl_leaf = leaf_by_name.get(cand.control)
    if conv_leaf is None or ctrl_leaf is None:
        return PairInfo(cand.convergent, cand.control)
    anc = tree.common_ancestor(conv_leaf, ctrl_leaf)

    convs: List[str] = []
    ctrls: List[str] = []
    for leaf in anc.get_terminals():
        nm = leaf.name or ""
        ph = pheno_for_algo.get(nm)
        if ph == 1:
            convs.append(nm)
        elif ph == -1:
            ctrls.append(nm)

    conv_choice = cand.convergent
    ctrl_choice = cand.control

    if method == "longest":
        if len(convs) > 1:
            conv_choice = _pick_longest(seq_lengths, convs)
        if len(ctrls) > 1:
            ctrl_choice = _pick_longest(seq_lengths, ctrls)
    elif method == "composite":
        duo = composite_duo.get((cand.convergent, cand.control))
        if duo:
            conv_choice, ctrl_choice = duo
    elif method == "contrast" and continuous:
        if len(convs) > 1:
            conv_choice = max(convs, key=lambda n: float(phenotypes.get(n, float("-inf"))))
        if len(ctrls) > 1:
            ctrl_choice = min(ctrls, key=lambda n: float(phenotypes.get(n, float("inf"))))
    elif method == "random":
        if len(convs) > 1:
            conv_choice = random.choice(convs)
        if len(ctrls) > 1:
            ctrl_choice = random.choice(ctrls)

    return PairInfo(conv_choice, ctrl_choice)


def auto_select_pairs(
    tree: Tree,
    phenotypes: Dict[str, float],
    *,
    method: str = "default",
    num_alternates: int = 0,
    max_combinations: int = 1,
    alignments_dir: str = "",
    lower_threshold: Optional[float] = None,
    upper_threshold: Optional[float] = None,
    quantile_tails_pct: Optional[float] = None,
    min_pct_diff: float = 0.0,
    seed: Optional[int] = None,
) -> List[PairInfo]:
    if seed is not None:
        random.seed(int(seed))

    leaf_names = {leaf.name or "" for leaf in tree.get_terminals()}
    phenotypes = {n: float(v) for n, v in phenotypes.items() if n in leaf_names}
    if not phenotypes:
        raise ValueError("No phenotype entries matched species present in the tree")

    method = (method or "default").strip().lower()
    if method == "simple_deterministic":
        method = "default"
    if method not in {"default", "longest", "shortest", "contrast", "composite", "random", "pct_contrast"}:
        raise ValueError(f"Unknown method: {method}")

    values = list(phenotypes.values())
    continuous = _pheno_is_continuous(values)

    parent_map = _build_parent_map(tree.root)
    leaf_by_name: Dict[str, Clade] = {}
    for leaf in tree.get_terminals():
        nm = leaf.name or ""
        if nm and nm not in leaf_by_name:
            leaf_by_name[nm] = leaf

    if method == "pct_contrast":
        if not continuous:
            raise ValueError("Method 'pct_contrast' is only available for continuous phenotypes")
        if any(float(v) <= 0.0 for v in values):
            raise ValueError("Method 'pct_contrast' requires strictly positive phenotype values")
        if int(num_alternates) > 0:
            raise ValueError("Alternates are not supported for method 'pct_contrast'")
        thr = max(0.0, float(min_pct_diff))
        candidates = _build_pct_contrast_candidates(tree, parent_map, leaf_by_name, phenotypes)
        selected_cands = _select_pct_contrast_candidates(candidates, parent_map, thr)
        if len(selected_cands) < 2:
            raise ValueError("ESL-PSC requires at least 2 valid contrast pairs")
        added = [PairInfo(c.convergent, c.control) for c in selected_cands]
        y_pos = _y_positions(tree)

        def anc_y(pair: PairInfo) -> float:
            conv_leaf = leaf_by_name.get(pair.convergent)
            ctrl_leaf = leaf_by_name.get(pair.control)
            if conv_leaf is None or ctrl_leaf is None:
                return 0.0
            try:
                anc = tree.common_ancestor(conv_leaf, ctrl_leaf)
            except Exception:
                return 0.0
            return float(y_pos.get(anc, 0.0))

        added.sort(key=anc_y)
        return added

    pheno_for_algo: Dict[str, int] = {}
    if continuous:
        if quantile_tails_pct is not None and float(quantile_tails_pct) > 0.0:
            lo, hi = thresholds_from_quantile_tails(values, float(quantile_tails_pct))
            lower_threshold, upper_threshold = float(lo), float(hi)
        if lower_threshold is None or upper_threshold is None:
            med = _median(values)
            lower_threshold, upper_threshold = float(med), float(med)
        low = float(lower_threshold)
        up = float(upper_threshold)
        if low > up:
            raise ValueError("Lower threshold must not exceed upper threshold")
        for name, val in phenotypes.items():
            if val > up:
                pheno_for_algo[name] = 1
            elif val < low:
                pheno_for_algo[name] = -1
    else:
        for name, val in phenotypes.items():
            if float(val) == 1.0:
                pheno_for_algo[name] = 1
            elif float(val) == -1.0:
                pheno_for_algo[name] = -1

    if sum(1 for v in pheno_for_algo.values() if v in (1, -1)) < 4:
        raise ValueError("Not enough species outside the thresholds for auto-selection")

    if num_alternates <= 0:
        max_combinations = 1

    if method in {"longest", "composite"} and not alignments_dir:
        raise ValueError("alignments_dir is required for method 'longest' and 'composite'")

    seq_lengths: Dict[str, int] = {}
    if method in {"longest", "composite"}:
        needed = [n for n, v in pheno_for_algo.items() if v in (1, -1)]
        seq_lengths = _scan_sequence_lengths(alignments_dir, needed)

    candidates: List[CandidatePair] = []
    prev_name: Optional[str] = None
    prev_ph = None
    prev_leaf: Optional[Clade] = None

    for leaf in list(tree.get_terminals()):
        name = leaf.name or ""
        ph = pheno_for_algo.get(name)
        if ph not in (1, -1):
            continue
        if prev_name is not None and ph != prev_ph:
            if prev_ph == 1:
                conv, ctrl = prev_name, name
                conv_leaf, ctrl_leaf = prev_leaf, leaf
            else:
                conv, ctrl = name, prev_name
                conv_leaf, ctrl_leaf = leaf, prev_leaf
            if conv_leaf is not None and ctrl_leaf is not None:
                anc = tree.common_ancestor(conv_leaf, ctrl_leaf)
                dist = _pair_distance(tree, parent_map, conv_leaf, ctrl_leaf)
                candidates.append(CandidatePair(conv, ctrl, dist, ancestor=anc))
        prev_name, prev_ph, prev_leaf = name, ph, leaf

    candidates.sort(key=lambda c: c.distance)

    composite_duo: Dict[Tuple[str, str], Tuple[str, str]] = {}
    composite_per_duo: Dict[Tuple[str, str], List[Tuple[str, str, float, float, float, float]]] = {}

    if method == "composite":
        trait_vals = [float(v) for v in phenotypes.values() if v is not None]
        med_trait = _median(trait_vals)
        devs = [abs(x - med_trait) for x in trait_vals]
        mad = _median(devs) if devs else 0.0
        p10 = _percentile(trait_vals, 0.10)
        p90 = _percentile(trait_vals, 0.90)
        central80 = p90 - p10
        floor = max(1e-12, 1e-9 * abs(med_trait))
        S_global = max(1.4826 * mad, (central80 / 2.563) if central80 > 0 else 0.0, floor)

        lens_all = list(seq_lengths.values())
        L_med = _median(lens_all) if lens_all else 0.0

        eps_abs = 0.15
        eps_rel = 0.05
        r_ok = 0.90
        r_bad = 0.50
        epsilon = 0.05

        for c in candidates:
            anc = c.ancestor
            if anc is None:
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

            per_duo: List[Tuple[str, str, float, float, float, float]] = []
            S_max = None
            for a in convs:
                for b in ctrls:
                    va = phenotypes.get(a)
                    vb = phenotypes.get(b)
                    if va is None or vb is None:
                        continue
                    diff = float(va) - float(vb)
                    if S_global <= 0.0:
                        continue
                    T = diff / S_global

                    la = seq_lengths.get(a)
                    lb = seq_lengths.get(b)
                    if la is None or lb is None or L_med <= 0.0:
                        G_len = epsilon
                        L_duo = 0.0
                    else:
                        try:
                            L_duo = 2.0 / (1.0 / float(la) + 1.0 / float(lb))
                        except Exception:
                            L_duo = 0.0
                        r = float(L_duo) / float(L_med) if L_med else 0.0
                        if r >= r_ok:
                            G_len = 1.0
                        elif r <= r_bad:
                            G_len = epsilon
                        else:
                            G_len = epsilon + (1.0 - epsilon) * ((r - r_bad) / (r_ok - r_bad))

                    S_val = T * G_len

                    a_leaf = leaf_by_name.get(a)
                    b_leaf = leaf_by_name.get(b)
                    dist = _duo_distance(tree, parent_map, a_leaf, b_leaf) if a_leaf and b_leaf else float("inf")

                    per_duo.append((a, b, S_val, dist, L_duo, T))
                    if S_max is None or S_val > S_max:
                        S_max = S_val

            if not per_duo or S_max is None:
                continue

            tie_threshold = max(S_max - eps_abs, (1.0 - eps_rel) * S_max)
            tie_set = [t for t in per_duo if t[2] >= tie_threshold]
            tie_set.sort(key=lambda t: (t[3], -t[5], -t[4], t[0], t[1]))
            a_best, b_best, _s, _d, _L, _T = tie_set[0]

            key = (c.convergent, c.control)
            composite_duo[key] = (a_best, b_best)
            composite_per_duo[key] = per_duo

    added: List[PairInfo] = []
    alt_scores: List[Tuple[int, str, str, float]] = []
    blocked: List[Clade] = []

    def overlaps_blocked(anc: Clade) -> bool:
        for b in blocked:
            if _is_descendant(parent_map, b, anc) or _is_descendant(parent_map, anc, b):
                return True
        return False

    for cand in candidates:
        if cand.ancestor is not None and overlaps_blocked(cand.ancestor):
            continue

        pair = _resolve_pair(
            tree,
            leaf_by_name,
            seq_lengths,
            cand,
            method,
            pheno_for_algo,
            phenotypes,
            composite_duo,
            continuous,
        )

        if num_alternates > 0:
            conv_leaf = leaf_by_name.get(pair.convergent)
            ctrl_leaf = leaf_by_name.get(pair.control)
            anc = None
            if conv_leaf is not None and ctrl_leaf is not None:
                try:
                    anc = tree.common_ancestor(conv_leaf, ctrl_leaf)
                except Exception:
                    anc = None
            if anc is not None:
                convs: List[str] = []
                ctrls: List[str] = []
                for leaf in anc.get_terminals():
                    nm = leaf.name or ""
                    ph = pheno_for_algo.get(nm)
                    if ph == 1 and nm != pair.convergent:
                        convs.append(nm)
                    elif ph == -1 and nm != pair.control:
                        ctrls.append(nm)

                def score_conv_alt(name: str) -> float:
                    other = pair.control
                    if method == "composite":
                        pd = composite_per_duo.get((cand.convergent, cand.control), [])
                        vals = [S for (a, b, S, _d, _L, _T) in pd if a == name and b == other]
                        return max(vals) if vals else float("-inf")
                    if method == "contrast" and continuous:
                        try:
                            return abs(float(phenotypes.get(name)) - float(phenotypes.get(other)))
                        except Exception:
                            return float("-inf")
                    if method == "longest":
                        la = float(seq_lengths.get(name, 0))
                        lb = float(seq_lengths.get(other, 0))
                        if la <= 0 or lb <= 0:
                            return 0.0
                        try:
                            return 2.0 / (1.0 / la + 1.0 / lb)
                        except Exception:
                            return 0.0
                    if method in ("shortest", "default"):
                        d = _pair_distance_names(tree, parent_map, leaf_by_name, name, other)
                        return -d
                    if method == "random":
                        return random.random()
                    d = _pair_distance_names(tree, parent_map, leaf_by_name, name, other)
                    return -d

                def score_ctrl_alt(name: str) -> float:
                    other = pair.convergent
                    if method == "composite":
                        pd = composite_per_duo.get((cand.convergent, cand.control), [])
                        vals = [S for (a, b, S, _d, _L, _T) in pd if a == other and b == name]
                        return max(vals) if vals else float("-inf")
                    if method == "contrast" and continuous:
                        try:
                            return abs(float(phenotypes.get(other)) - float(phenotypes.get(name)))
                        except Exception:
                            return float("-inf")
                    if method == "longest":
                        la = float(seq_lengths.get(other, 0))
                        lb = float(seq_lengths.get(name, 0))
                        if la <= 0 or lb <= 0:
                            return 0.0
                        try:
                            return 2.0 / (1.0 / la + 1.0 / lb)
                        except Exception:
                            return 0.0
                    if method in ("shortest", "default"):
                        d = _pair_distance_names(tree, parent_map, leaf_by_name, other, name)
                        return -d
                    if method == "random":
                        return random.random()
                    d = _pair_distance_names(tree, parent_map, leaf_by_name, other, name)
                    return -d

                if convs:
                    convs_scored = [(n, score_conv_alt(n)) for n in convs]
                    convs_scored = [t for t in convs_scored if t[1] != float("-inf")]
                    convs_scored.sort(key=lambda t: t[1], reverse=True)
                    chosen = [n for (n, _s) in convs_scored[: int(num_alternates)]]
                    pair.conv_alts.extend(chosen)
                    idx_in_added = len(added)
                    for n, s in convs_scored[: int(num_alternates)]:
                        alt_scores.append((idx_in_added, "conv", n, float(s)))

                if ctrls:
                    ctrls_scored = [(n, score_ctrl_alt(n)) for n in ctrls]
                    ctrls_scored = [t for t in ctrls_scored if t[1] != float("-inf")]
                    ctrls_scored.sort(key=lambda t: t[1], reverse=True)
                    chosen = [n for (n, _s) in ctrls_scored[: int(num_alternates)]]
                    pair.ctrl_alts.extend(chosen)
                    idx_in_added = len(added)
                    for n, s in ctrls_scored[: int(num_alternates)]:
                        alt_scores.append((idx_in_added, "ctrl", n, float(s)))

        added.append(pair)
        if cand.ancestor is not None:
            blocked.append(cand.ancestor)

    if len(added) < 2:
        raise ValueError("ESL-PSC requires at least 2 valid contrast pairs")

    if max_combinations is not None and int(max_combinations) > 0 and alt_scores:
        tmp_pairs = [PairInfo(p.convergent, p.control, list(p.conv_alts), list(p.ctrl_alts)) for p in added]
        alt_scores_sorted = sorted(alt_scores, key=lambda x: (x[3], x[1], x[2]))
        while _compute_combo_count(tmp_pairs, cap=int(max_combinations)) > int(max_combinations) and alt_scores_sorted:
            idx, role, name, _score = alt_scores_sorted.pop(0)
            if idx < 0 or idx >= len(tmp_pairs):
                continue
            if role == "conv" and name in tmp_pairs[idx].conv_alts:
                tmp_pairs[idx].conv_alts.remove(name)
            elif role == "ctrl" and name in tmp_pairs[idx].ctrl_alts:
                tmp_pairs[idx].ctrl_alts.remove(name)
        for i, p in enumerate(tmp_pairs):
            added[i].conv_alts = p.conv_alts
            added[i].ctrl_alts = p.ctrl_alts

    ancestors: List[Optional[Clade]] = []
    for p in added:
        conv_leaf = leaf_by_name.get(p.convergent)
        ctrl_leaf = leaf_by_name.get(p.control)
        if conv_leaf is None or ctrl_leaf is None:
            ancestors.append(None)
            continue
        try:
            ancestors.append(tree.common_ancestor(conv_leaf, ctrl_leaf))
        except Exception:
            ancestors.append(None)

    keep: List[PairInfo] = []
    for i, anc_i in enumerate(ancestors):
        if anc_i is None:
            continue
        nested = False
        for j, anc_j in enumerate(ancestors):
            if i == j or anc_j is None:
                continue
            if _is_descendant(parent_map, anc_i, anc_j):
                nested = True
                break
        if not nested:
            keep.append(added[i])
    added = keep

    y_pos = _y_positions(tree)

    def anc_y(pair: PairInfo) -> float:
        conv_leaf = leaf_by_name.get(pair.convergent)
        ctrl_leaf = leaf_by_name.get(pair.control)
        if conv_leaf is None or ctrl_leaf is None:
            return 0.0
        try:
            anc = tree.common_ancestor(conv_leaf, ctrl_leaf)
        except Exception:
            return 0.0
        return float(y_pos.get(anc, 0.0))

    added.sort(key=anc_y)

    return added


def write_species_groups(pairs: List[PairInfo], output_path: str) -> None:
    if not pairs:
        raise ValueError("No pairs to write")
    with open(output_path, "w", encoding="utf-8") as f:
        for pair in pairs:
            conv_line = ",".join([pair.convergent] + list(pair.conv_alts))
            ctrl_line = ",".join([pair.control] + list(pair.ctrl_alts))
            f.write(f"{conv_line}\n{ctrl_line}\n")
