"""Fast scanning of alignments for CCS and control convergence."""
from __future__ import annotations

import os
from collections import Counter
from typing import Callable, Dict, List

from gui.core.fasta_io import read_fasta


def _parse_species_groups(path: str) -> List[tuple[List[str], List[str]]]:
    """Return list of (convergent_species, control_species) combos."""
    combos: List[tuple[List[str], List[str]]] = []
    if not path or not os.path.exists(path):
        return combos
    lines: List[str] = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    for i in range(0, len(lines), 2):
        conv_line = lines[i]
        ctrl_line = lines[i + 1] if i + 1 < len(lines) else ""
        conv = [sp.strip() for sp in conv_line.replace(",", " ").split() if sp.strip()]
        ctrl = [sp.strip() for sp in ctrl_line.replace(",", " ").split() if sp.strip()]
        combos.append((conv, ctrl))
    return combos


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


def fast_scan_alignments(
    alignment_dir: str,
    species_groups_file: str,
    outgroup_species: str,
    progress_cb: Callable[[int, int], None] | None = None,
) -> List[Dict[str, float]]:
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
    results: List[Dict[str, float]] = []
    total = len(files)
    for idx, fname in enumerate(files, 1):
        path = os.path.join(alignment_dir, fname)
        records = read_fasta(path)
        species_seq = {sp: seq for sp, seq in records}
        if not species_seq:
            if progress_cb:
                progress_cb(idx, total)
            continue
        seq_len = len(next(iter(species_seq.values())))
        true_counts: List[int] = []
        ctrl_counts: List[int] = []
        for conv, ctrl in combos:
            ccs = 0
            ctrl_conv = 0
            for pos in range(seq_len):
                conv_aa = [species_seq.get(sp, "?")[pos] for sp in conv]
                ctrl_aa = [species_seq.get(sp, "?")[pos] for sp in ctrl]
                out_aa = [species_seq.get(outgroup_species, "?")[pos]]
                clean_conv = [x for x in conv_aa if x not in ("?", "-")]
                clean_ctrl = [x for x in ctrl_aa if x not in ("?", "-")]
                clean_out = [x for x in out_aa if x not in ("?", "-")]
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
                if clean_out and len(set(clean_out)) == 1:
                    out_res = clean_out[0]
                    if clean_conv and all(r == out_res for r in clean_conv):
                        ctrl_counter = Counter(clean_ctrl)
                        for res, cnt in ctrl_counter.items():
                            if res != out_res and cnt >= 2:
                                ctrl_conv += 1
                                break
            true_counts.append(ccs)
            ctrl_counts.append(ctrl_conv)
        avg_true = sum(true_counts) / len(combos) if combos else 0.0
        avg_ctrl = sum(ctrl_counts) / len(combos) if combos else 0.0
        results.append(
            {
                "gene": os.path.splitext(fname)[0],
                "avg_true": avg_true,
                "avg_control": avg_ctrl,
                "diff": avg_true - avg_ctrl,
            }
        )
        if progress_cb:
            progress_cb(idx, total)
    results.sort(key=lambda x: x["avg_true"], reverse=True)
    return results
