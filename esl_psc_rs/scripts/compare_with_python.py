#!/usr/bin/env python3
import argparse
import csv
import math
import os
import subprocess
import tempfile
from pathlib import Path


def pearson(xs, ys):
    if not xs:
        return float("nan")
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return 0.0
    return num / (dx * dy)


def load_gene_gss(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[row["gene_name"]] = float(row["highest_gss"])
    return out


def load_sites(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[(row["gene_name"], int(row["position"]))] = float(row["pss"])
    return out


def main():
    parser = argparse.ArgumentParser(description="Compare unified Rust output vs Python ESL-PSC baseline")
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument(
        "--response-matrix",
        default="esl_psc_output/photo_single_LC_matrix_species_groups_response_matrices/combo_0.txt",
    )
    parser.add_argument("--alignments-dir", default="demo_data/photosynthesis/alignments")
    parser.add_argument("--l1-min", type=float, default=0.1)
    parser.add_argument("--l1-max", type=float, default=0.2)
    parser.add_argument("--l2-min", type=float, default=0.1)
    parser.add_argument("--l2-max", type=float, default=0.2)
    parser.add_argument("--lambda-step", type=float, default=0.1)
    parser.add_argument("--penalty-type", default="std")
    parser.add_argument("--min-gene-pearson", type=float, default=0.9)
    parser.add_argument("--min-site-pearson", type=float, default=0.9)
    parser.add_argument("--min-top10-overlap", type=int, default=5)
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    response = (root / args.response_matrix).resolve()
    alignments = (root / args.alignments_dir).resolve()

    if not response.exists():
        raise SystemExit(f"response matrix not found: {response}")
    if not alignments.exists():
        raise SystemExit(f"alignments dir not found: {alignments}")

    with tempfile.TemporaryDirectory(prefix="esl_py_") as py_out, tempfile.TemporaryDirectory(prefix="esl_rs_") as rs_out:
        py_out = Path(py_out)
        rs_out = Path(rs_out)

        py_cmd = [
            "python3",
            "-m",
            "esl_psc_cli.esl_integrator",
            "--response_matrix_path",
            str(response),
            "--input_alignments_dir",
            str(alignments),
            "--preprocessed_dir_name",
            "py_cmp_combo0",
            "--output_file_base_name",
            "pybase",
            "--output_dir",
            str(py_out),
            "--esl_inputs_outputs_dir",
            str(py_out / "pre"),
            "--initial_lambda1",
            str(args.l1_min),
            "--final_lambda1",
            str(args.l1_max),
            "--initial_lambda2",
            str(args.l2_min),
            "--final_lambda2",
            str(args.l2_max),
            "--lambda_step",
            str(args.lambda_step),
            "--group_penalty_type",
            args.penalty_type,
            "--no_pred_output",
            "--show_selected_sites",
        ]

        subprocess.run(py_cmd, cwd=root, check=True)
        # Python baseline creates alignments/paths.txt as a side effect.
        # Remove it to keep the repository clean after comparison runs.
        maybe_paths = alignments / "paths.txt"
        if maybe_paths.exists():
            maybe_paths.unlink()

        rs_cmd = [
            "cargo",
            "run",
            "--quiet",
            "--release",
            "--",
            "--response-matrix-path",
            str(response),
            "--input-alignments-dir",
            str(alignments),
            "--output-file-base-name",
            "rustbase",
            "--output-dir",
            str(rs_out),
            "--initial-lambda1",
            str(args.l1_min),
            "--final-lambda1",
            str(args.l1_max),
            "--initial-lambda2",
            str(args.l2_min),
            "--final-lambda2",
            str(args.l2_max),
            "--lambda-step",
            str(args.lambda_step),
            "--group-penalty-type",
            args.penalty_type,
            "--no-pred-output",
            "--show-selected-sites",
        ]

        subprocess.run(rs_cmd, cwd=root / "esl_psc_rs", check=True)

        py_gene = load_gene_gss(py_out / "pybase_gene_ranks.csv")
        rs_gene = load_gene_gss(rs_out / "rustbase_gene_ranks.csv")
        common_genes = sorted(set(py_gene) & set(rs_gene))
        gx = [py_gene[g] for g in common_genes]
        gy = [rs_gene[g] for g in common_genes]
        gene_pear = pearson(gx, gy)

        py_sites = load_sites(py_out / "pybase_selected_sites.csv")
        rs_sites = load_sites(rs_out / "rustbase_selected_sites.csv")
        common_sites = sorted(set(py_sites) & set(rs_sites))
        sx = [py_sites[k] for k in common_sites]
        sy = [rs_sites[k] for k in common_sites]
        site_pear = pearson(sx, sy)

        py_top10 = [k for k, _ in sorted(py_gene.items(), key=lambda kv: kv[1], reverse=True)[:10]]
        rs_top10 = [k for k, _ in sorted(rs_gene.items(), key=lambda kv: kv[1], reverse=True)[:10]]
        top10_overlap = len(set(py_top10) & set(rs_top10))

        print(f"gene_pearson={gene_pear:.6f} (common_genes={len(common_genes)})")
        print(f"selected_site_pearson={site_pear:.6f} (common_sites={len(common_sites)})")
        print(f"top10_overlap={top10_overlap}")

        ok = (
            gene_pear >= args.min_gene_pearson
            and site_pear >= args.min_site_pearson
            and top10_overlap >= args.min_top10_overlap
        )
        if not ok:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
