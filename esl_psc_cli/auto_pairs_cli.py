from __future__ import annotations

import argparse
import os
import sys


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

    args = parser.parse_args(argv)

    tree_file = os.path.abspath(args.tree_file)
    pheno_file = os.path.abspath(args.species_pheno_path)
    output_path = os.path.abspath(args.output_path)

    if not os.path.isfile(tree_file):
        print(f"Error: tree file not found: {tree_file}", file=sys.stderr)
        return 2
    if not os.path.isfile(pheno_file):
        print(f"Error: phenotype file not found: {pheno_file}", file=sys.stderr)
        return 2

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    try:
        tree = load_tree(tree_file)
        phenos = load_phenotypes(pheno_file)
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
            seed=args.seed,
        )
        write_species_groups(pairs, output_path)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
