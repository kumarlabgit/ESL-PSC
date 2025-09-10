"""Fast scanning of alignments for CCS, control convergence, and CS-derived significance."""
from __future__ import annotations

import os
import math
from collections import Counter
from typing import Callable, Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

from gui.core.fasta_io import read_fasta
from esl_psc_cli.deletion_canceler import (
    parse_species_groups as cli_parse_species_groups,
)


def _parse_species_groups(path: str) -> List[tuple[List[str], List[str]]]:
    """Return list of (convergent_species, control_species) combos derived
    from the CLI's species-groups parser to ensure exact consistency.

    This wraps ``esl_psc_cli.deletion_canceler.parse_species_groups`` (which
    returns tuples of species chosen across all lines) and splits each tuple
    by index parity into Convergent (even indices) and Control (odd indices).
    """
    try:
        raw_combos = cli_parse_species_groups(path)  # list[tuple[str, ...]]
    except Exception:
        return []
    conv_ctrl: List[tuple[List[str], List[str]]] = []
    for tup in raw_combos:
        picks = list(tup)
        conv = [picks[i] for i in range(0, len(picks), 2)]
        ctrl = [picks[i] for i in range(1, len(picks), 2)]
        conv_ctrl.append((conv, ctrl))
    return conv_ctrl


def list_species(alignment_dir: str) -> List[str]:
    """Collect all species names present across alignment files."""
    species: set[str] = set()
    if not alignment_dir or not os.path.isdir(alignment_dir):
        return []
    for fname in os.listdir(alignment_dir):
        if not fname.endswith(".fas"):
            continue
        for rec_id, _ in read_fasta(os.path.join(alignment_dir, fname)):
            species.add(rec_id)
    return sorted(species)


# ─────────────────────────────────────────────────────────────────────────────
# Per-file worker for parallel fast scan
# ─────────────────────────────────────────────────────────────────────────────
def _scan_file_worker(
    alignment_dir: str,
    fname: str,
    combos: List[Tuple[List[str], List[str]]],
    outgroup_species: str,
) -> Dict[str, float | int]:
    path = os.path.join(alignment_dir, fname)
    records = read_fasta(path)
    species_seq = {sp: seq for sp, seq in records}
    if not species_seq:
        return {
            "gene": os.path.splitext(fname)[0],
            "avg_true": 0.0,
            "avg_control": 0.0,
            "diff": 0.0,
            "variable_sites": 0,
            "cs_sites_ge_4": 0,
            "k_pairs": 0,
        }

    # Use the maximum sequence length seen to be robust to minor length differences
    seq_len = max((len(s) for s in species_seq.values()), default=0)

    # Safe accessor: return '?' when species is missing or position is out of range
    def get_aa(species: str, pos: int) -> str:
        seq = species_seq.get(species)
        if not seq or pos < 0 or pos >= len(seq):
            return "?"
        try:
            return seq[pos]
        except Exception:
            return "?"

    # No CS p-value helpers – we only compute CS≥2 site counts for speed

    files_true_counts: List[int] = []
    files_ctrl_counts: List[int] = []
    true_den = 0
    ctrl_den = 0

    # Track max CS≥2 across combos (report the combo with the most CS sites)
    best_cs_sites = 0
    best_k_pairs = 0

    # Pre-compute total variable sites across ALL species in this file
    # Mirror esl_psc_functions.is_variable_site logic for consistency
    variable_sites_total = 0
    species_list_all = list(species_seq.keys())
    for pos in range(seq_len):
        aa_all = [get_aa(sp, pos) for sp in species_list_all]
        counted = Counter(aa_all)
        if '-' in counted:
            del counted['-']
        num_left = sum(counted.values())
        if len(counted) == 0 or len(counted) == 1:
            continue
        if len(counted) == num_left:
            continue  # all singletons
        if max(counted.values()) == num_left - 1:
            continue  # one singleton
        variable_sites_total += 1

    # Track per-combo true convergence and diff (true - control) for cross-combo ranking later
    per_combo_true: List[float | None] = [None] * len(combos)
    per_combo_diff: List[float | None] = [None] * len(combos)

    for combo_idx, (conv_group, ctrl_group) in enumerate(combos):
        # Build explicit pairs by index and keep only pairs with both species present
        indexed_pairs = list(zip(conv_group, ctrl_group))
        pairs_present = [(a, b) for (a, b) in indexed_pairs if a in species_seq and b in species_seq]
        conv_present = [a for (a, _b) in pairs_present]
        ctrl_present = [b for (_a, b) in pairs_present]
        eligible_true = len(conv_present) >= 2 and len(ctrl_present) >= 1
        eligible_ctrl = len(ctrl_present) >= 2 and len(conv_present) >= 1
        if not (eligible_true or eligible_ctrl):
            continue

        ccs = 0
        ctrl_conv = 0
        cs_sites = 0

        # Pre-extract sequence strings for this combo to avoid dict lookups per position
        conv_seqs = [species_seq[a] for (a, _b) in pairs_present]
        ctrl_seqs = [species_seq[b] for (_a, b) in pairs_present]
        out_seq = species_seq.get(outgroup_species, "")

        for pos in range(seq_len):
            # Raw residues at this position
            conv_aa = [s[pos] if 0 <= pos < len(s) else '?' for s in conv_seqs]
            ctrl_aa = [s[pos] if 0 <= pos < len(s) else '?' for s in ctrl_seqs]
            out_aa = [out_seq[pos]] if 0 <= pos < len(out_seq) else ['?']

            # Mask singletons among Convergent+Control
            cc_counts = Counter(conv_aa + ctrl_aa)
            conv_ns = conv_aa[:]
            ctrl_ns = ctrl_aa[:]
            for lst in (conv_ns, ctrl_ns):
                for i, r in enumerate(lst):
                    if r != '-' and cc_counts.get(r, 0) == 1:
                        lst[i] = '?'

            # Clean lists for checks
            clean_conv = [x for x in conv_ns if x not in ('?', '-')]
            clean_ctrl = [x for x in ctrl_ns if x not in ('?', '-')]
            clean_out = [x for x in out_aa if x not in ('?', '-')]

            # CCS detection (only meaningful if eligible_true)
            if eligible_true:
                if (
                    clean_ctrl
                    and clean_out
                    and len(set(clean_ctrl)) == 1
                    and len(set(clean_out)) == 1
                    and list(set(clean_ctrl))[0] == list(set(clean_out))[0]
                ):
                    ctrl_res = clean_ctrl[0]
                    conv_counter = Counter(clean_conv)
                    for res, cnt in conv_counter.items():
                        if res != ctrl_res and cnt >= 2:
                            ccs += 1
                            break

            # Control-convergence detection (only meaningful if eligible_ctrl)
            if eligible_ctrl and clean_out and len(set(clean_out)) == 1:
                out_res = clean_out[0]
                if clean_conv and all(r == out_res for r in clean_conv):
                    ctrl_counter = Counter(clean_ctrl)
                    for res, cnt in ctrl_counter.items():
                        if res != out_res and cnt >= 2:
                            ctrl_conv += 1
                            break

            # CS computation per site (final_score only)
            diff_counts = Counter(clean_conv)
            diff_counts.subtract(clean_ctrl)
            raw_score = sum(abs(v) for v in diff_counts.values())
            gap_count = sum(1 for a in conv_ns if a == '-') + sum(1 for b in ctrl_ns if b == '-')
            final_score = raw_score - gap_count
            if final_score < 0:
                final_score = 0
            if final_score >= 4:
                cs_sites += 1

        if eligible_true:
            files_true_counts.append(ccs)
            true_den += 1
            per_combo_true[combo_idx] = float(ccs)
        if eligible_ctrl:
            files_ctrl_counts.append(ctrl_conv)
            ctrl_den += 1

        # Per-combo diff only if both components are eligible
        if eligible_true and eligible_ctrl:
            per_combo_diff[combo_idx] = float(ccs - ctrl_conv)

        # Track the combo with the most CS≥2 sites
        if cs_sites > best_cs_sites:
            best_cs_sites = cs_sites
            best_k_pairs = len(pairs_present)

    avg_true = (sum(files_true_counts) / true_den) if true_den else 0.0
    avg_ctrl = (sum(files_ctrl_counts) / ctrl_den) if ctrl_den else 0.0

    return {
        "gene": os.path.splitext(fname)[0],
        "avg_true": avg_true,
        "avg_control": avg_ctrl,
        "diff": avg_true - avg_ctrl,
        "variable_sites": variable_sites_total,
        "cs_sites_ge_4": int(best_cs_sites),
        "k_pairs": int(best_k_pairs),
        "per_combo_true": per_combo_true,
        "per_combo_diff": per_combo_diff,
    }


def fast_scan_alignments(
    alignment_dir: str,
    species_groups_file: str,
    outgroup_species: str,
    progress_cb: Callable[[int, int], None] | None = None,
    response_dir: str | None = None,
    n_jobs: int | None = None,
    top_frac: float = 0.01,
) -> List[Dict[str, float | int]]:
    """Scan all alignments and compute convergence metrics per gene.

    Parameters
    ----------
    alignment_dir: str
        Directory containing ``.fas`` alignment files.
    species_groups_file: str
        Path to species groups definition.
    outgroup_species: str
        Species to treat as the outgroup.
    progress_cb: callable, optional
        Callback receiving ``(current, total)`` file counts.
    """
    combos = _parse_species_groups(species_groups_file)
    files = sorted([f for f in os.listdir(alignment_dir) if f.endswith(".fas")])
    total = len(files)
    results: List[Dict[str, float | int]] = []

    # Default n_jobs if not specified
    if n_jobs is None:
        try:
            cpu = os.cpu_count() or 1
            n_jobs = min(4, cpu)
        except Exception:
            n_jobs = 1

    # Serial path
    if n_jobs <= 1:
        for idx, fname in enumerate(files, 1):
            res = _scan_file_worker(alignment_dir, fname, combos, outgroup_species)
            results.append(res)
            if progress_cb:
                progress_cb(idx, total)
        # Post-process combo-based ranking if multiple combos
        if len(combos) > 1:
            _apply_combo_top_ranking(results, len(combos), top_frac)
            _apply_combo_top_ranking_by_diff(results, len(combos), top_frac)
            results.sort(
                key=lambda x: (
                    x.get("num_combos_top_frac_by_diff", 0),
                    x.get("num_combos_top_frac", 0),
                    x.get("avg_true", 0.0),
                ),
                reverse=True,
            )
        else:
            results.sort(key=lambda x: x["avg_true"], reverse=True)
        return results

    # Parallel path using processes
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        future_map = {ex.submit(_scan_file_worker, alignment_dir, fname, combos, outgroup_species): fname for fname in files}
        done = 0
        for fut in as_completed(future_map):
            try:
                res = fut.result()
                results.append(res)
            except Exception:
                # On failure for a file, skip but keep going
                pass
            done += 1
            if progress_cb:
                progress_cb(done, total)
    if len(combos) > 1:
        _apply_combo_top_ranking(results, len(combos), top_frac)
        _apply_combo_top_ranking_by_diff(results, len(combos), top_frac)
        results.sort(
            key=lambda x: (
                x.get("num_combos_top_frac_by_diff", 0),
                x.get("num_combos_top_frac", 0),
                x.get("avg_true", 0.0),
            ),
            reverse=True,
        )
    else:
        results.sort(key=lambda x: x["avg_true"], reverse=True)
    return results

def _apply_combo_top_ranking(results: List[Dict[str, float | int]], n_combos: int, top_frac: float) -> None:
    """For multiple combos, compute per-gene count of combos where the gene is in the
    top fraction by true convergence for that combo. Mutates results to add:
      - num_combos_top_frac
      - top_fraction (float)
    """
    # Build per-combo lists of (index, value)
    per_combo_lists: List[List[tuple[int, float]]] = [[] for _ in range(n_combos)]
    for i, row in enumerate(results):
        vals = row.get("per_combo_true") or []
        # Ensure vals length matches n_combos
        if len(vals) < n_combos:
            vals = list(vals) + [None] * (n_combos - len(vals))
        for j, v in enumerate(vals):
            if v is None:
                continue
            try:
                per_combo_lists[j].append((i, float(v)))
            except Exception:
                continue

    in_top_counts = [0] * len(results)
    for j in range(n_combos):
        lst = per_combo_lists[j]
        if not lst:
            continue
        # Sort descending by value
        lst.sort(key=lambda t: t[1], reverse=True)
        # Determine top-k by fraction, include ties at boundary
        k = max(1, int(len(lst) * top_frac))
        if k >= len(lst):
            chosen = set(idx for idx, _ in lst)
        else:
            cutoff_val = lst[k - 1][1]
            chosen = set(idx for idx, val in lst if val >= cutoff_val)
        for idx in chosen:
            in_top_counts[idx] += 1

    for i, row in enumerate(results):
        row["num_combos_top_frac"] = in_top_counts[i]
        row["top_fraction"] = float(top_frac)

def _apply_combo_top_ranking_by_diff(results: List[Dict[str, float | int]], n_combos: int, top_frac: float) -> None:
    """Like _apply_combo_top_ranking but uses per-combo diffs (true - control).
    Adds:
      - num_combos_top_frac_by_diff
      - top_fraction_by_diff
    """
    per_combo_lists: List[List[tuple[int, float]]] = [[] for _ in range(n_combos)]
    for i, row in enumerate(results):
        vals = row.get("per_combo_diff") or []
        if len(vals) < n_combos:
            vals = list(vals) + [None] * (n_combos - len(vals))
        for j, v in enumerate(vals):
            if v is None:
                continue
            try:
                per_combo_lists[j].append((i, float(v)))
            except Exception:
                continue
    in_top_counts = [0] * len(results)
    for j in range(n_combos):
        lst = per_combo_lists[j]
        if not lst:
            continue
        lst.sort(key=lambda t: t[1], reverse=True)
        k = max(1, int(len(lst) * top_frac))
        if k >= len(lst):
            chosen = set(idx for idx, _ in lst)
        else:
            cutoff_val = lst[k - 1][1]
            chosen = set(idx for idx, val in lst if val >= cutoff_val)
        for idx in chosen:
            in_top_counts[idx] += 1
    for i, row in enumerate(results):
        row["num_combos_top_frac_by_diff"] = in_top_counts[i]
        row["top_fraction_by_diff"] = float(top_frac)
