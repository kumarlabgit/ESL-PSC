# esl-psc_cli

`esl-psc_cli` is a unified Rust executable that runs:

1. Multimatrix orchestration (`species_groups_file` or `response_dir`)
2. Gap-cancel + preprocess in memory
3. Sparse-group-lasso model fitting over lambda grids (in memory)
4. ESL-style outputs (gene ranks, selected sites, species predictions, model files)

This removes the repeated Python/C file handoff loop (`paths.txt`, preprocess dirs,
group-penalty rewrites, repeated solver input/output parsing) and keeps each combo
in memory while sweeping lambdas/penalties.

## Scope

Supported run modes:

- Single response matrix (`--response-matrix-path`)
- Multi response matrix directory (`--response-dir`)
- Species-groups multimatrix (`--species-groups-file`)

Important options:

- Alignment paths: `--alignments-dir` / `--input-alignments-dir`, `--prediction-alignments-dir`
- Hyperparameters: lambda ranges, grid/logspace, group penalty settings
- Gap-cancel controls: `--use-uncanceled-alignments`, `--cancel-only-partner`,
  `--cancel-tri-allelic`, `--min-pairs`, `--outgroup-species`
- Outputs: `--no-pred-output`, `--no-genes-output`, `--show-selected-sites`

## Build

```bash
cd esl_unified_rs
cargo build --release
```

## Run

Example (multimatrix):

```bash
./target/release/esl-psc_cli \
  --species-groups-file ../photo_multi_species_groups.txt \
  --alignments-dir ../photosynthesis_alignments \
  --prediction-alignments-dir ../photosynthesis_alignments \
  --species-pheno-path ../photo_species_phenotypes.txt \
  --output-file-base-name rust_mm \
  --output-dir /tmp/esl_unified_mm \
  --initial-lambda1 0.1 --final-lambda1 0.2 \
  --initial-lambda2 0.1 --final-lambda2 0.2 \
  --lambda-step 0.1 \
  --group-penalty-type std \
  --show-selected-sites
```

Outputs:

- `<base>_gene_ranks.csv`
- `<base>_selected_sites.csv` (if enabled)
- `<base>_species_predictions.csv` (unless `--no-pred-output`)
- `models_unified_rs/*.txt` (non-zero feature weights + intercept)

## Baseline comparison

Single-response parity check against Python `esl_integrator`:

```bash
python3 scripts/compare_with_python.py
```

Multimatrix parity + timing benchmark against Python `esl_multimatrix`:

```bash
python3 scripts/benchmark_multimatrix_vs_python.py
```

The multimatrix benchmark reports:

- Pearson correlation of multimatrix gene scores (`highest_ever_gss`)
- Pearson correlation of selected-site PSS on common sites
- Prediction Pearson correlation on aligned prediction rows
- Top-10 overlap
- Wall-clock speedup (Python / Rust)

and exits non-zero if thresholds are below thresholds.

## CLI Surface Parity Check

You can verify flag-surface parity against Python ESL-PSC with:

```bash
python3 scripts/check_cli_parity.py
```

This checks that all Python CLI long options are present in `esl-psc_cli`.
