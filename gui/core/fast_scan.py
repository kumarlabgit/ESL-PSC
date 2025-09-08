"""Fast scanning of alignments for CCS and control convergence."""
from __future__ import annotations

import os
from collections import Counter
from typing import Callable, Dict, List

from gui.core.fasta_io import read_fasta
from esl_psc_cli.deletion_canceler import (
    parse_species_groups as cli_parse_species_groups,
)
from esl_psc_cli.esl_psc_functions import count_var_sites
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def fast_scan_alignments(
    alignment_dir: str,
    species_groups_file: str,
    outgroup_species: str,
    progress_cb: Callable[[int, int], None] | None = None,
    response_dir: str | None = None,
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
    results: List[Dict[str, float]] = []
    total = len(files)
    for idx, fname in enumerate(files, 1):
        path = os.path.join(alignment_dir, fname)
        records = read_fasta(path)
        seq_records = [SeqRecord(Seq(seq), id=sp) for sp, seq in records]
        species_seq = {sp: seq for sp, seq in records}
        if not species_seq:
            if progress_cb:
                progress_cb(idx, total)
            continue
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
        # Compute CCS/Control-convergence per combo, then average across eligible combos
        true_counts: List[int] = []
        ctrl_counts: List[int] = []
        true_den = 0  # combos with >=2 convergent and >=1 control present
        ctrl_den = 0  # combos with >=2 control and >=1 convergent present
        for combo_index, (conv_group, ctrl_group) in enumerate(combos):
            conv_present = [sp for sp in conv_group if sp in species_seq]
            ctrl_present = [sp for sp in ctrl_group if sp in species_seq]
            # Determine eligibility for each metric separately
            eligible_true = len(conv_present) >= 2 and len(ctrl_present) >= 1
            eligible_ctrl = len(ctrl_present) >= 2 and len(conv_present) >= 1
            if not (eligible_true or eligible_ctrl):
                continue
            ccs = 0
            ctrl_conv = 0
            for pos in range(seq_len):
                # Raw residues at this position
                conv_aa = [get_aa(sp, pos) for sp in conv_present]
                ctrl_aa = [get_aa(sp, pos) for sp in ctrl_present]
                out_aa = [get_aa(outgroup_species, pos)]

                # Mask singletons among Convergent+Control (match SiteViewer logic)
                # Count on raw residues in Convergent+Control (including '-')
                cc_counts = Counter(conv_aa + ctrl_aa)
                conv_ns = conv_aa[:]
                ctrl_ns = ctrl_aa[:]
                for lst in (conv_ns, ctrl_ns):
                    for i, r in enumerate(lst):
                        if r != '-' and cc_counts.get(r, 0) == 1:
                            lst[i] = '?'

                # Clean up for detection
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
            if eligible_true:
                true_counts.append(ccs)
                true_den += 1
            if eligible_ctrl:
                ctrl_counts.append(ctrl_conv)
                ctrl_den += 1
            

        avg_true = (sum(true_counts) / true_den) if true_den else 0.0
        avg_ctrl = (sum(ctrl_counts) / ctrl_den) if ctrl_den else 0.0
        var_sites = count_var_sites(seq_records)
        results.append(
            {
                "gene": os.path.splitext(fname)[0],
                "avg_true": avg_true,
                "avg_control": avg_ctrl,
                "diff": avg_true - avg_ctrl,
                "variable_sites": var_sites,
            }
        )
        if progress_cb:
            progress_cb(idx, total)
    results.sort(key=lambda x: x["avg_true"], reverse=True)
    return results
