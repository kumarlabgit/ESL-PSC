"""Fast scanning of alignments for CCS, control convergence, and CS-derived significance.

Optionally uses a Rust backend binary (fast_scan_rs) if available to accelerate
per-file scanning. Falls back to the Python implementation otherwise.
"""
from __future__ import annotations

import os
import json
import subprocess
import sys
from collections import Counter
from typing import Callable, Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

from gui.core.fasta_io import read_fasta
from esl_psc_cli.deletion_canceler import (
    parse_species_groups as cli_parse_species_groups,
)
from esl_psc_cli import esl_psc_functions as ecf


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


def _parse_species_groups_two_pair(path: str) -> List[tuple[List[str], List[str]]]:
    """Return all 2x2 combos implied by a species groups file.

    Interpret the groups file as alternating convergent/control lines per pair.
    For N pairs, generate combos for every choice of two distinct pairs (i, j),
    i < j. If a line contains multiple comma-separated species, generate a
    variant for each option.

    Example for 3 pairs (lines C1, K1, C2, K2, C3, K3):
      combos include (C1,K1,C2,K2), (C1,K1,C3,K3), (C2,K2,C3,K3) with expansion
      across any multi-species alternatives per line.
    """
    try:
        with open(path, encoding="utf-8", errors="ignore") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
    except Exception:
        return []
    if not lines or len(lines) % 2 != 0:
        return []
    groups: List[List[str]] = []
    for ln in lines:
        opts = [sp.strip() for sp in ln.split(",") if sp.strip()]
        if not opts:
            return []
        groups.append(opts)
    n_pairs = len(groups) // 2
    combos: List[tuple[List[str], List[str]]] = []
    # Choose 2 pairs (i, j)
    for i in range(n_pairs):
        for j in range(i + 1, n_pairs):
            conv_i = groups[2 * i]
            ctrl_i = groups[2 * i + 1]
            conv_j = groups[2 * j]
            ctrl_j = groups[2 * j + 1]
            for a in conv_i:
                for b in ctrl_i:
                    for c in conv_j:
                        for d in ctrl_j:
                            combos.append(([a, c], [b, d]))
    return combos


def list_species(alignment_dir: str, max_no_new: int = 200) -> List[str]:
    """Collect all species names present across alignment files.
    
    Optimized to scan only FASTA header lines (">...") to avoid reading full
    sequences, which improves responsiveness when choosing an outgroup.
    
    Early-stop heuristic: if no new species are discovered across
    ``max_no_new`` consecutive files (default 200), assume the set is complete
    and stop scanning. Pass ``max_no_new=0`` to disable early stopping.
    """
    species: set[str] = set()
    if not alignment_dir or not os.path.isdir(alignment_dir):
        return []
    files = sorted([f for f in os.listdir(alignment_dir) if ecf.is_fasta(f)])
    no_new_streak = 0
    for fname in files:
        fpath = os.path.join(alignment_dir, fname)
        before = len(species)
        try:
            with open(fpath, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line or line[0] != ">":
                        continue
                    # Take the first token after '>' as species ID (consistent with read_fasta)
                    header = line[1:].strip()
                    if not header:
                        continue
                    sp = header.split()[0]
                    if sp:
                        species.add(sp)
        except Exception:
            # Skip unreadable or malformed files
            pass
        after = len(species)
        if after == before:
            no_new_streak += 1
        else:
            no_new_streak = 0
        if max_no_new and no_new_streak >= max_no_new:
            break
    return sorted(species)


# ─────────────────────────────────────────────────────────────────────────────
# Per-file worker for parallel fast scan
# ─────────────────────────────────────────────────────────────────────────────
def _scan_file_worker(
    alignment_dir: str,
    fname: str,
    combos: List[Tuple[List[str], List[str]]],
    outgroup_species: str,
    min_out_ctrl_agreement: float,
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
                # Require non-gap, non-missing outgroup residue and at least the
                # configured fraction of control residues to match the outgroup.
                if clean_out:
                    out_uniform = list(set(clean_out))
                    if len(out_uniform) == 1:
                        out_res = out_uniform[0]
                        if clean_ctrl:
                            matches = sum(1 for r in clean_ctrl if r == out_res)
                            total = len(clean_ctrl)
                            agree = (matches / total) if total > 0 else 0.0
                            if agree >= max(0.0, min(1.0, float(min_out_ctrl_agreement))):
                                # Count CCS if any convergent residue != out_res occurs at least twice
                                conv_counter = Counter(clean_conv)
                                for res, cnt in conv_counter.items():
                                    if res != out_res and cnt >= 2:
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
    two_pair_combos: bool = False,
    min_out_ctrl_agreement: float = 1.0,
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
    combos = _parse_species_groups_two_pair(species_groups_file) if two_pair_combos else _parse_species_groups(species_groups_file)
    files = sorted([f for f in os.listdir(alignment_dir) if ecf.is_fasta(f)])
    total = len(files)
    results: List[Dict[str, float | int]] = []

    # Default n_jobs if not specified
    if n_jobs is None:
        try:
            cpu = os.cpu_count() or 1
            n_jobs = min(4, cpu)
        except Exception:
            n_jobs = 1

    # Try Rust backend if available
    rs_bin = _detect_fast_scan_rs()
    # If fractional control-outgroup agreement is requested, prefer a freshly built
    # target/release binary (newer) over any packaged bin/ version. If not found,
    # fall back to Python to ensure correct semantics.
    try:
        if float(min_out_ctrl_agreement) != 1.0:
            repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            cand = os.path.join(repo_root, "fast_scan_rs", "target", "release", "fast_scan_rs")
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                rs_bin = cand
            else:
                rs_bin = None
    except Exception:
        rs_bin = None
    if rs_bin:
        try:
            # Single Rust invocation; if a progress callback is provided, we stream
            # progress from the Rust process stderr in real time.
            rs_results = _run_fast_scan_rs(
                rs_bin,
                alignment_dir,
                files,
                combos,
                outgroup_species,
                cs_threshold=4,
                min_out_ctrl_agreement=float(min_out_ctrl_agreement),
                progress_cb=progress_cb,
                total=total,
                done_offset=0,
            )
            results = rs_results
            if len(combos) > 1:
                _apply_combo_top_ranking(results, len(combos), top_frac)
                _apply_combo_top_ranking_by_diff(results, len(combos), top_frac)
                _apply_combo_top_ranking_by_ratio(results, len(combos), top_frac)
                results.sort(
                    key=lambda x: (
                        x.get("num_combos_top_frac", 0),
                        x.get("num_combos_top_frac_by_ratio", 0),
                        x.get("num_combos_top_frac_by_diff", 0),
                        x.get("avg_true", 0.0),
                        x.get("gene", ""),
                    ),
                    reverse=True,
                )
            else:
                results.sort(key=lambda x: x["avg_true"], reverse=True)
            return results
        except Exception:
            # Fallback to Python path on any error
            pass

    # Serial or parallel path (Python)
    if n_jobs <= 1:
        for idx, fname in enumerate(files, 1):
            res = _scan_file_worker(alignment_dir, fname, combos, outgroup_species, float(min_out_ctrl_agreement))
            results.append(res)
            if progress_cb:
                progress_cb(idx, total)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            future_map = {ex.submit(_scan_file_worker, alignment_dir, fname, combos, outgroup_species, float(min_out_ctrl_agreement)): fname for fname in files}
            done = 0
            for fut in as_completed(future_map):
                try:
                    res = fut.result()
                    results.append(res)
                except Exception:
                    # Skip failed file but continue
                    pass
                done += 1
                if progress_cb:
                    progress_cb(done, total)
    # Post-process combo-based ranking if multiple combos
    if len(combos) > 1:
        _apply_combo_top_ranking(results, len(combos), top_frac)
        _apply_combo_top_ranking_by_diff(results, len(combos), top_frac)
        _apply_combo_top_ranking_by_ratio(results, len(combos), top_frac)
        results.sort(
            key=lambda x: (
                x.get("num_combos_top_frac", 0),
                x.get("num_combos_top_frac_by_ratio", 0),
                x.get("num_combos_top_frac_by_diff", 0),
                x.get("avg_true", 0.0),
                x.get("gene", ""),
            ),
            reverse=True,
        )
    else:
        results.sort(key=lambda x: x["avg_true"], reverse=True)
    return results


def _detect_fast_scan_rs() -> str | None:
    """Return path to Rust fast_scan_rs binary if available, else None.

    Detection priority:
      1) Environment variable FAST_SCAN_RS pointing to the binary
      2) bin/fast_scan_rs (relative to repo root)
      3) fast_scan_rs/target/release/fast_scan_rs (relative to repo root)
      4) bin/ next to the executable (packaged app/onefile)
    """
    # Allow disabling via env for testing or debugging
    if os.environ.get("FAST_SCAN_RS_DISABLE", "0") in {"1", "true", "True"}:
        return None
    # 1) Env var
    env_path = os.environ.get("FAST_SCAN_RS")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path
    # Compute repo root from this file (editable/source check)
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    # 2) bin/fast_scan_rs[_mac] relative to repo root (editable install)
    if sys.platform == "darwin":
        cand = os.path.join(repo_root, "bin", "fast_scan_rs_mac")
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    cand = os.path.join(repo_root, "bin", "fast_scan_rs")
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    # 3) fast_scan_rs/target/release/fast_scan_rs (relative to repo root)
    cand = os.path.join(repo_root, "fast_scan_rs", "target", "release", "fast_scan_rs")
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    # Also try Windows .exe (even on WSL)
    cand = os.path.join(repo_root, "bin", "fast_scan_rs.exe")
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    # 4) Packaged app path: bin/ next to the executable
    try:
        exe_dir = os.path.abspath(os.path.dirname(getattr(sys, 'executable', sys.argv[0])))
        # Prefer mac-suffixed name on macOS bundles
        if sys.platform == "darwin":
            cand = os.path.join(exe_dir, "bin", "fast_scan_rs_mac")
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return cand
        cand = os.path.join(exe_dir, "bin", "fast_scan_rs")
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
        cand = os.path.join(exe_dir, "bin", "fast_scan_rs.exe")
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    except Exception:
        pass
    return None


def _run_fast_scan_rs(
    bin_path: str,
    alignment_dir: str,
    files: List[str],
    combos: List[Tuple[List[str], List[str]]],
    outgroup_species: str,
    cs_threshold: int = 4,
    min_out_ctrl_agreement: float | None = None,
    progress_cb: Callable[[int, int], None] | None = None,
    total: int | None = None,
    done_offset: int = 0,
) -> List[Dict[str, float | int]]:
    """Invoke the Rust fast_scan_rs binary and parse its JSON output.

    Returns a list of rows with the same schema as the Python worker.
    """
    # Build JSON spec
    spec = {
        "alignment_dir": alignment_dir,
        "files": files,
        "combos": [{"conv": c, "ctrl": d} for (c, d) in combos],
        "outgroup": outgroup_species,
        "cs_threshold": int(cs_threshold),
        "emit_progress": bool(progress_cb is not None and total),
    }
    # Backward compatibility: only include the field if it differs from default 1.0
    try:
        if min_out_ctrl_agreement is not None and float(min_out_ctrl_agreement) != 1.0:
            spec["min_out_ctrl_agreement"] = float(min_out_ctrl_agreement)
    except Exception:
        pass
    # Use a temporary file for stdout to avoid pipe backpressure while we read stderr for progress
    import tempfile
    with tempfile.TemporaryFile(mode="w+b") as tmp_out:
        try:
            proc = subprocess.Popen(
                [bin_path],
                stdin=subprocess.PIPE,
                stdout=tmp_out,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
        except Exception as e:
            raise RuntimeError(f"fast_scan_rs failed to start: {e}")
        # Send JSON spec and close stdin
        try:
            assert proc.stdin is not None
            proc.stdin.write(json.dumps(spec))
            proc.stdin.close()
        except Exception:
            pass
        # Stream stderr for progress lines like: PROGRESS <n>
        done = 0
        if proc.stderr is not None and progress_cb and total:
            for line in proc.stderr:
                if not line:
                    break
                line = line.strip()
                if line.startswith("PROGRESS "):
                    try:
                        n = int(line.split(" ")[1])
                        done = n
                        progress_cb(min(done_offset + done, total), total)
                    except Exception:
                        pass
        # Wait for completion
        rc = proc.wait()
        # Read stdout from temp file
        tmp_out.seek(0)
        stdout = tmp_out.read().decode("utf-8", errors="replace").strip()
        if rc != 0:
            raise RuntimeError(f"fast_scan_rs exited with code {rc}. Stderr may contain details.")
    if not stdout:
        raise RuntimeError("fast_scan_rs produced no output")
    try:
        data = json.loads(stdout)
    except Exception as e:
        raise RuntimeError(f"fast_scan_rs invalid JSON: {e}\nFirst 200 chars: {stdout[:200]}")
    # Validate and coerce to expected list of dicts
    if not isinstance(data, list):
        raise RuntimeError("fast_scan_rs output was not a list")
    out: List[Dict[str, float | int]] = []
    for row in data:
        if not isinstance(row, dict):
            continue
        # Ensure required keys exist; default if missing
        out.append({
            "gene": row.get("gene", ""),
            "avg_true": float(row.get("avg_true", 0.0)),
            "avg_control": float(row.get("avg_control", 0.0)),
            "diff": float(row.get("diff", 0.0)),
            "variable_sites": int(row.get("variable_sites", 0) or 0),
            "cs_sites_ge_4": int(row.get("cs_sites_ge_4", 0) or 0),
            "k_pairs": int(row.get("k_pairs", 0) or 0),
            "per_combo_true": row.get("per_combo_true", []),
            "per_combo_diff": row.get("per_combo_diff", []),
        })
    return out

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

def _apply_combo_top_ranking_by_ratio(results: List[Dict[str, float | int]], n_combos: int, top_frac: float) -> None:
    """Compute per-gene count of combos where the true/(control+1) ratio is in the top fraction.

    Uses per-combo true counts and diffs (true - control) to derive control.
    If per-combo diff is missing, assume control=0 for that combo when true is present.
    Adds:
      - num_combos_top_frac_by_ratio
      - top_fraction_by_ratio
    """
    per_combo_lists: List[List[tuple[int, float]]] = [[] for _ in range(n_combos)]
    for i, row in enumerate(results):
        tvals = row.get("per_combo_true") or []
        dvals = row.get("per_combo_diff") or []
        # Ensure lengths match n_combos, padding with None
        if len(tvals) < n_combos:
            tvals = list(tvals) + [None] * (n_combos - len(tvals))
        if len(dvals) < n_combos:
            dvals = list(dvals) + [None] * (n_combos - len(dvals))
        for j in range(n_combos):
            tv = tvals[j]
            if tv is None:
                continue
            try:
                t = float(tv)
            except Exception:
                continue
            dv = dvals[j]
            ctrl = 0.0
            if dv is not None:
                try:
                    d = float(dv)
                    ctrl = max(0.0, t - d)
                except Exception:
                    ctrl = 0.0
            ratio = t / (ctrl + 1.0)
            per_combo_lists[j].append((i, ratio))

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
        row["num_combos_top_frac_by_ratio"] = in_top_counts[i]
        row["top_fraction_by_ratio"] = float(top_frac)
