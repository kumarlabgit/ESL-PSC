from __future__ import annotations

import argparse
import os
import sys


def _build_output_paths(output_path: str, num_random_sets: int, output_is_dir: bool = False) -> list[str]:
    if num_random_sets <= 1:
        return [output_path]

    # If a directory path is requested, write numbered files inside it.
    if output_is_dir:
        output_dir = output_path
        stem, ext = "auto_pairs_groups", ".txt"
    else:
        output_dir = os.path.dirname(output_path)
        base_name = os.path.basename(output_path)
        stem, ext = os.path.splitext(base_name)
        if not stem:
            stem = "auto_pairs_groups"

    width = max(3, len(str(int(num_random_sets))))
    return [
        os.path.join(output_dir, f"{stem}_{i:0{width}d}{ext}")
        for i in range(1, int(num_random_sets) + 1)
    ]


def main(argv=None) -> int:
    try:
        from esl_psc_cli.auto_pairs import (
            auto_select_pairs,
            load_phenotypes,
            load_tree,
            write_species_groups,
        )
    except ModuleNotFoundError as exc:
        if getattr(exc, "name", "") in {"Bio", "Bio.Phylo"}:
            print(
                "Error: biopython is required for auto pair selection. Install it (e.g. `pip install biopython`).",
                file=sys.stderr,
            )
            return 1
        raise

    parser = argparse.ArgumentParser(
        description="Generate an ESL-PSC species groups file via automatic contrast pair selection."
    )
    parser.add_argument(
        "--tree_file",
        required=True,
        help="Newick or NEXUS tree file",
    )
    parser.add_argument(
        "--species_pheno_path",
        required=True,
        help="Species phenotype file (CSV: species,value) supporting binary or continuous values",
    )
    parser.add_argument(
        "--output_path",
        required=True,
        help="Output species groups file path",
    )
    parser.add_argument(
        "--method",
        default="default",
        choices=["default", "longest", "shortest", "contrast", "composite", "random"],
        help="Tie-breaking method for ambiguous ancestor clades",
    )
    parser.add_argument(
        "--num_alternates",
        type=int,
        default=0,
        help="Number of alternates to include per side of each pair (default: 0)",
    )
    parser.add_argument(
        "--max_combinations",
        type=int,
        default=1,
        help="Maximum number of total combinations across all pairs (default: 1)",
    )
    parser.add_argument(
        "--alignments_dir",
        default="",
        help="Alignments directory (required for --method longest or composite)",
    )
    parser.add_argument(
        "--lower_threshold",
        type=float,
        default=None,
        help="Lower threshold for continuous phenotypes (values below => control)",
    )
    parser.add_argument(
        "--upper_threshold",
        type=float,
        default=None,
        help="Upper threshold for continuous phenotypes (values above => convergent)",
    )
    parser.add_argument(
        "--quantile_tails_pct",
        type=float,
        default=None,
        help="If set (>0), override thresholds using symmetric quantile tails in percent (0-50)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (used only for --method random)",
    )
    parser.add_argument(
        "--num_random_sets",
        type=int,
        default=1,
        help=(
            "Generate N random species-groups files with numbered names in the "
            "output location (only with --method random)"
        ),
    )

    args = parser.parse_args(argv)

    tree_file = os.path.abspath(args.tree_file)
    pheno_file = os.path.abspath(args.species_pheno_path)
    output_arg = str(args.output_path)
    output_path = os.path.abspath(output_arg)
    output_is_dir = (
        output_arg.endswith(os.sep)
        or (os.altsep is not None and output_arg.endswith(os.altsep))
    )

    if not os.path.isfile(tree_file):
        print(f"Error: tree file not found: {tree_file}", file=sys.stderr)
        return 2
    if not os.path.isfile(pheno_file):
        print(f"Error: phenotype file not found: {pheno_file}", file=sys.stderr)
        return 2

    if int(args.num_random_sets) < 1:
        parser.error("--num_random_sets must be >= 1")
    if int(args.num_random_sets) > 1 and str(args.method) != "random":
        parser.error("--num_random_sets > 1 requires --method random")

    output_paths = _build_output_paths(
        output_path,
        int(args.num_random_sets),
        output_is_dir=bool(output_is_dir),
    )
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    for pth in output_paths:
        pth_dir = os.path.dirname(pth)
        if pth_dir:
            os.makedirs(pth_dir, exist_ok=True)

    try:
        tree = load_tree(tree_file)
        phenos = load_phenotypes(pheno_file)
        for i, out_path in enumerate(output_paths):
            run_seed = None if args.seed is None else int(args.seed) + i
            pairs = auto_select_pairs(
                tree,
                phenos,
                method=str(args.method),
                num_alternates=int(args.num_alternates),
                max_combinations=int(args.max_combinations),
                alignments_dir=str(args.alignments_dir or ""),
                lower_threshold=args.lower_threshold,
                upper_threshold=args.upper_threshold,
                quantile_tails_pct=args.quantile_tails_pct,
                seed=run_seed,
            )
            write_species_groups(pairs, out_path)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
