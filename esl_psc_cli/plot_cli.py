from __future__ import annotations

import argparse
import os
import sys

from esl_psc_cli import esl_psc_functions as ecf


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate ESL-PSC prediction plots from a species predictions CSV."
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["kde", "violin", "continuous"],
        help="Plot mode: kde, violin, or continuous",
    )
    parser.add_argument(
        "--pred_csv",
        required=True,
        help="Path to <output_base>_species_predictions.csv",
    )
    parser.add_argument(
        "--title",
        required=True,
        help="Plot title/output base name",
    )
    parser.add_argument(
        "--min_genes",
        type=int,
        default=0,
        help="Minimum genes threshold (default: 0)",
    )
    parser.add_argument(
        "--pheno_name1",
        default=None,
        help="Positive phenotype display name (binary plots only)",
    )
    parser.add_argument(
        "--pheno_name2",
        default=None,
        help="Negative phenotype display name (binary plots only)",
    )
    args = parser.parse_args(argv)

    pred_csv = os.path.abspath(args.pred_csv)
    if not os.path.isfile(pred_csv):
        print(f"Error: predictions CSV not found: {pred_csv}", file=sys.stderr)
        return 2

    try:
        if args.mode == "continuous":
            ecf.continuous_pred_plot(
                pred_csv,
                args.title,
                min_genes=int(args.min_genes),
            )
        else:
            pheno_names = None
            if args.pheno_name1 and args.pheno_name2:
                pheno_names = (str(args.pheno_name1), str(args.pheno_name2))
            ecf.rmse_range_pred_plots(
                pred_csv,
                args.title,
                pheno_names=pheno_names,
                min_genes=int(args.min_genes),
                plot_type=str(args.mode),
            )
    except Exception as exc:
        print(f"Error: plot generation failed: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
