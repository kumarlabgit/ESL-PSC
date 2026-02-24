import argparse
import json
import os
import sys
from typing import List, Dict

# Reuse the GUI core fast_scan implementation for parsing, Rust detection, and Python fallback.
from gui.core import fast_scan as fs
from esl_psc_cli import esl_psc_functions as ecf


def _postprocess_and_sort(results: List[Dict], n_combos: int, top_frac: float) -> List[Dict]:
    """Mirror gui.core.fast_scan.fast_scan_alignments postprocessing and sorting."""
    out = list(results)
    if n_combos > 1:
        fs._apply_combo_top_ranking(out, n_combos, top_frac)
        fs._apply_combo_top_ranking_by_diff(out, n_combos, top_frac)
        fs._apply_combo_top_ranking_by_ratio(out, n_combos, top_frac)
        out.sort(
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
        out.sort(key=lambda x: x.get("avg_true", 0.0), reverse=True)
    return out


def _print_progress(cur: int, total: int) -> None:
    # Simple CLI progress: write a single line that updates in place
    try:
        pct = (100.0 * cur / total) if total else 0.0
        sys.stderr.write(f"\rScanning: {cur}/{total} ({pct:5.1f}%)")
        sys.stderr.flush()
        if cur >= total:
            sys.stderr.write("\n")
            sys.stderr.flush()
    except Exception:
        pass


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Site Counter alignments for CCS and related metrics (prefers Rust backend)."
    )
    parser.add_argument(
        "--alignments_dir",
        required=True,
        help="Directory containing FASTA alignment files (.fas, .fasta, .fa, .faa)",
    )
    parser.add_argument(
        "--species_groups_file",
        required=True,
        help="Path to species groups file (alternating convergent/control lines)",
    )
    parser.add_argument(
        "--outgroup_species",
        required=False,
        help="Outgroup species identifier (required unless using --tree_file)",
    )
    parser.add_argument(
        "--tree_file",
        type=str,
        help="Newick or NEXUS tree file for parsimony ancestral reconstruction",
    )
    parser.add_argument(
        "--require_unambiguous_mrca",
        action="store_true",
        help="Require unambiguous ancestral residues (exclude ambiguous sites from CCS detection)",
    )
    parser.add_argument(
        "--compute_mrca_representative",
        action="store_true",
        help="Compute a single representative ancestral residue sequence (disabled by default; not needed for set-based CCS logic).",
    )
    parser.add_argument(
        "--two_pair_combos",
        action="store_true",
        help="Interpret species groups by generating all 2x2 pair combos",
    )
    parser.add_argument(
        "--min_out_ctrl_agreement",
        type=float,
        default=1.0,
        help="Minimum fraction [0,1] of control residues matching the outgroup to count CCS (default: 1.0)",
    )
    parser.add_argument(
        "--top_frac",
        type=float,
        default=0.01,
        help="Top fraction per combo used for cross-combo ranking (default: 0.01)",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Output CSV file path",
    )

    args = parser.parse_args(argv)

    align_dir = os.path.abspath(args.alignments_dir)
    groups_path = os.path.abspath(args.species_groups_file)
    outgroup = args.outgroup_species
    tree_file = args.tree_file
    two_pair = bool(args.two_pair_combos)
    require_unamb = bool(args.require_unambiguous_mrca)
    compute_rep = bool(args.compute_mrca_representative)
    min_agree = float(args.min_out_ctrl_agreement)
    top_frac = float(args.top_frac)
    
    # Validate that either outgroup or tree is provided
    if not outgroup and not tree_file:
        print("Error: must provide either --outgroup_species or --tree_file", file=sys.stderr)
        return 2
    if tree_file:
        tree_file = os.path.abspath(tree_file)
        if not os.path.isfile(tree_file):
            print(f"Error: tree file not found: {tree_file}", file=sys.stderr)
            return 2
        # Set placeholder outgroup for ancestral mode
        if not outgroup:
            outgroup = "ANCESTRAL_MRCA"

    # Validate inputs
    if not os.path.isdir(align_dir):
        print(f"Error: alignments directory not found: {align_dir}", file=sys.stderr)
        return 2
    if not os.path.isfile(groups_path):
        print(f"Error: species groups file not found: {groups_path}", file=sys.stderr)
        return 2

    # Build file list and combos once
    files = sorted([f for f in os.listdir(align_dir) if ecf.is_fasta(f)])
    if not files:
        print("Error: no FASTA files (.fas, .fasta, .fa, .faa) found in alignments directory", file=sys.stderr)
        return 2

    combos = (
        fs._parse_species_groups_two_pair(groups_path)
        if two_pair
        else fs._parse_species_groups(groups_path)
    )
    if not combos:
        print("Error: no valid species combos could be constructed from the groups file", file=sys.stderr)
        return 2
    n_combos = len(combos)
    
    # Collect analysis species for ancestral reconstruction
    analysis_species = None
    if tree_file:
        analysis_species = set()
        for conv, ctrl in combos:
            analysis_species.update(conv)
            analysis_species.update(ctrl)

    results: List[Dict]
    used_rust = False

    # Attempt Rust path unless disabled
    rs_bin = fs._detect_site_counter_rs()
    # If fractional agreement or tree is requested, prefer a newer target build if present
    if rs_bin and (min_agree != 1.0 or tree_file):
        try:
            repo_root = os.path.abspath(os.path.join(os.path.dirname(fs.__file__), "..", ".."))
            cand = os.path.join(repo_root, "fast_scan_rs", "target", "release", "site_counter_rs")
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                rs_bin = cand
            else:
                rs_bin = None
        except Exception:
            rs_bin = None

    if rs_bin:
        try:
            results = fs._run_site_counter_rs(
                rs_bin,
                align_dir,
                files,
                combos,
                outgroup,
                cs_threshold=4,
                min_out_ctrl_agreement=min_agree,
                progress_cb=_print_progress,
                total=len(files),
                done_offset=0,
                tree_file=tree_file,
                analysis_species=analysis_species,
                require_unambiguous_mrca=require_unamb,
                compute_mrca_representative=compute_rep,
            )
            used_rust = True
            results = _postprocess_and_sort(results, n_combos, top_frac)
        except Exception as e:
            # Warn and fall back to Python path
            print(
                f"Warning: Rust Site Counter backend failed ({e}). Falling back to Python implementation...",
                file=sys.stderr,
            )
            # Ensure Python fallback is used for this invocation
            os.environ["SITE_COUNTER_RS_DISABLE"] = "1"
            results = fs.fast_scan_alignments(
                align_dir,
                groups_path,
                outgroup,
                progress_cb=_print_progress,
                two_pair_combos=two_pair,
                min_out_ctrl_agreement=min_agree,
                tree_file=tree_file,
                require_unambiguous_mrca=require_unamb,
                compute_mrca_representative=compute_rep,
            )
    else:
        # Direct Python path (or forced)
        os.environ["SITE_COUNTER_RS_DISABLE"] = "1"
        results = fs.fast_scan_alignments(
            align_dir,
            groups_path,
            outgroup,
            progress_cb=_print_progress,
            two_pair_combos=two_pair,
            min_out_ctrl_agreement=min_agree,
            tree_file=tree_file,
            require_unambiguous_mrca=require_unamb,
            compute_mrca_representative=compute_rep,
        )

    # Emit CSV output (required path)
    keys = []
    seen = set()
    for row in results:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    lines = [",".join(keys)]
    def _fmt(v):
        if isinstance(v, float):
            return f"{v:.6g}"
        if isinstance(v, (list, tuple)):
            try:
                return json.dumps(v)
            except Exception:
                return str(v)
        return str(v)
    for row in results:
        lines.append(
            ",".join(_fmt(row.get(k, "")) for k in keys)
        )
    text = "\n".join(lines)
    with open(args.output_path, "w", encoding="utf-8") as fh:
        fh.write(text + ("\n" if not text.endswith("\n") else ""))
    # Print a small footer noting the backend used (stderr so it doesn't corrupt stdout)
    sys.stderr.write(f"Finished Site Counter using {'Rust' if used_rust else 'Python'} backend.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
