#!/usr/bin/env python3
import argparse
import csv
import math
import shutil
import subprocess
import tempfile
import time
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


def norm_float(s):
    return round(float(s), 6)


def load_gene_metric(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[row["gene_name"]] = float(row["highest_ever_gss"])
    return out


def load_sites(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[(row["gene_name"], int(row["position"]))] = float(row["pss"])
    return out


def load_predictions(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (
                row["species_combo"],
                norm_float(row["lambda1"]),
                norm_float(row["lambda2"]),
                norm_float(row["penalty_term"]),
                row["species"],
            )
            out[key] = float(row["SPS"])
    return out


def timed_run(cmd, cwd):
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.perf_counter() - t0


def main():
    p = argparse.ArgumentParser(description="Benchmark + parity check: Python esl_multimatrix vs unified Rust")
    p.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    p.add_argument("--species-groups-file", default="photo_multi_species_groups.txt")
    p.add_argument("--alignments-dir", default="photosynthesis_alignments")
    p.add_argument("--prediction-alignments-dir", default="photosynthesis_alignments")
    p.add_argument("--species-pheno-path", default="photo_species_phenotypes.txt")
    p.add_argument("--l1-min", type=float, default=0.1)
    p.add_argument("--l1-max", type=float, default=0.2)
    p.add_argument("--l2-min", type=float, default=0.1)
    p.add_argument("--l2-max", type=float, default=0.2)
    p.add_argument("--lambda-step", type=float, default=0.1)
    p.add_argument("--penalty-type", default="std")
    p.add_argument("--repeats", type=int, default=1)
    p.add_argument("--skip-predictions", action="store_true")
    p.add_argument("--keep-temp", action="store_true")

    p.add_argument("--min-gene-pearson", type=float, default=0.9)
    p.add_argument("--min-site-pearson", type=float, default=0.9)
    p.add_argument("--min-pred-pearson", type=float, default=0.9)
    p.add_argument("--min-top10-overlap", type=int, default=5)
    p.add_argument("--min-speedup", type=float, default=5.0)

    args = p.parse_args()

    root = Path(args.repo_root).resolve()
    species_groups_file = (root / args.species_groups_file).resolve()
    alignments_dir = (root / args.alignments_dir).resolve()
    pred_alignments_dir = (root / args.prediction_alignments_dir).resolve()
    species_pheno_path = (root / args.species_pheno_path).resolve()

    if not species_groups_file.exists():
        raise SystemExit(f"species groups file not found: {species_groups_file}")

    temp_parent = Path(tempfile.mkdtemp(prefix="esl_mm_bench_"))
    py_out = temp_parent / "python_out"
    rs_out = temp_parent / "rust_out"
    py_out.mkdir(parents=True, exist_ok=True)
    rs_out.mkdir(parents=True, exist_ok=True)

    py_cmd = [
        "python3",
        "-m",
        "esl_psc_cli.esl_multimatrix",
        "--species_groups_file",
        str(species_groups_file),
        "--alignments_dir",
        str(alignments_dir),
        "--prediction_alignments_dir",
        str(pred_alignments_dir),
        "--species_pheno_path",
        str(species_pheno_path),
        "--output_dir",
        str(py_out),
        "--output_file_base_name",
        "py_mm",
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
        "--show_selected_sites",
        "--no_checkpoint",
    ]

    rust_bin = root / "esl_unified_rs" / "target" / "release" / "esl-psc_cli"
    if not rust_bin.exists():
        rust_bin = root / "esl_unified_rs" / "target" / "release" / "esl_unified_rs"
    rs_cmd = [
        str(rust_bin),
        "--species-groups-file",
        str(species_groups_file),
        "--alignments-dir",
        str(alignments_dir),
        "--prediction-alignments-dir",
        str(pred_alignments_dir),
        "--species-pheno-path",
        str(species_pheno_path),
        "--output-dir",
        str(rs_out),
        "--output-file-base-name",
        "rs_mm",
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
        "--show-selected-sites",
    ]

    if args.skip_predictions:
        py_cmd.append("--no_pred_output")
        rs_cmd.append("--no-pred-output")

    py_times = []
    rs_times = []

    py_times.append(timed_run(py_cmd, root))
    rs_times.append(timed_run(rs_cmd, root))

    for rep in range(1, args.repeats):
        py_rep = temp_parent / f"python_out_{rep}"
        rs_rep = temp_parent / f"rust_out_{rep}"
        py_rep.mkdir(parents=True, exist_ok=True)
        rs_rep.mkdir(parents=True, exist_ok=True)

        py_rep_cmd = py_cmd.copy()
        py_rep_cmd[py_rep_cmd.index("--output_dir") + 1] = str(py_rep)
        rs_rep_cmd = rs_cmd.copy()
        rs_rep_cmd[rs_rep_cmd.index("--output-dir") + 1] = str(rs_rep)

        py_times.append(timed_run(py_rep_cmd, root))
        rs_times.append(timed_run(rs_rep_cmd, root))

    py_mean = sum(py_times) / len(py_times)
    rs_mean = sum(rs_times) / len(rs_times)
    speedup = py_mean / rs_mean if rs_mean > 0 else float("inf")

    py_gene = load_gene_metric(py_out / "py_mm_gene_ranks.csv")
    rs_gene = load_gene_metric(rs_out / "rs_mm_gene_ranks.csv")
    common_genes = sorted(set(py_gene) & set(rs_gene))
    gene_pear = pearson([py_gene[g] for g in common_genes], [rs_gene[g] for g in common_genes])

    py_sites = load_sites(py_out / "py_mm_selected_sites.csv")
    rs_sites = load_sites(rs_out / "rs_mm_selected_sites.csv")
    common_sites = sorted(set(py_sites) & set(rs_sites))
    site_pear = pearson([py_sites[k] for k in common_sites], [rs_sites[k] for k in common_sites])

    py_top10 = [k for k, _ in sorted(py_gene.items(), key=lambda kv: kv[1], reverse=True)[:10]]
    rs_top10 = [k for k, _ in sorted(rs_gene.items(), key=lambda kv: kv[1], reverse=True)[:10]]
    top10_overlap = len(set(py_top10) & set(rs_top10))

    print(f"python_times={py_times}")
    print(f"rust_times={rs_times}")
    print(f"python_mean={py_mean:.6f}s")
    print(f"rust_mean={rs_mean:.6f}s")
    print(f"speedup={speedup:.2f}x")

    print(f"gene_pearson={gene_pear:.6f} (common_genes={len(common_genes)})")
    print(f"selected_site_pearson={site_pear:.6f} (common_sites={len(common_sites)})")
    print(f"top10_overlap={top10_overlap}")
    print(f"site_counts_python={len(py_sites)} site_counts_rust={len(rs_sites)}")

    if not args.skip_predictions:
        py_pred = load_predictions(py_out / "py_mm_species_predictions.csv")
        rs_pred = load_predictions(rs_out / "rs_mm_species_predictions.csv")
        common_pred = sorted(set(py_pred) & set(rs_pred))
        pred_pear = pearson([py_pred[k] for k in common_pred], [rs_pred[k] for k in common_pred])
        print(f"prediction_pearson={pred_pear:.6f} (common_predictions={len(common_pred)})")
    else:
        pred_pear = float("nan")

    ok = (
        gene_pear >= args.min_gene_pearson
        and site_pear >= args.min_site_pearson
        and top10_overlap >= args.min_top10_overlap
        and speedup >= args.min_speedup
    )

    if not args.skip_predictions:
        ok = ok and pred_pear >= args.min_pred_pearson

    if not args.keep_temp:
        shutil.rmtree(temp_parent, ignore_errors=True)
    else:
        print(f"kept_temp_dir={temp_parent}")

    if not ok:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
