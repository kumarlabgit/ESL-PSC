# Auto Pair Selection CLI

This CLI generates an ESL-PSC species groups file from a phylogenetic tree and a species phenotype file.

Recommended command:

```bash
esl-psc pairs --help
```

Python module form (equivalent):

```bash
python -m esl_psc_cli.auto_pairs_cli --help
```

## Inputs

- `--tree_file`: Newick or NEXUS tree file.
- `--species_pheno_path`: CSV with `species,value` rows.
  - Binary values (`-1`, `1`) are supported.
  - Continuous values are supported.
- `--output_path`: Path for the output species groups file.

Output format:
- 2 lines per contrast pair.
- Odd line = convergent side.
- Even line = control side.
- If alternates are requested, they are comma-delimited on the same line.

## Basic Usage

```bash
esl-psc pairs \
  --tree_file /path/tree.nwk \
  --species_pheno_path /path/phenotypes.csv \
  --output_path /path/species_groups.txt
```

## Options

| Flag | Type | Default | Notes |
|---|---|---:|---|
| `--method` | string | `simple_deterministic` | One of: `simple_deterministic`, `default` (alias), `longest`, `shortest`, `contrast`, `composite`, `random`, `pct_contrast`. |
| `--num_alternates` | int | `0` | Number of alternates per side of each pair. |
| `--max_combinations` | int | `1` | Maximum total combinations across all pairs. |
| `--alignments_dir` | path | `""` | Required for `--method longest` or `--method composite`. |
| `--lower_threshold` | float | auto | Continuous mode: values below are controls (`-1`). |
| `--upper_threshold` | float | auto | Continuous mode: values above are convergent (`+1`). |
| `--quantile_tails_pct` | float | none | Continuous mode: overrides thresholds using symmetric quantile tails (0-50). |
| `--min_pct_diff` | float | `0.0` | Used by `--method pct_contrast` (positive continuous phenotypes only). |
| `--seed` | int | none | Random seed (used with `--method random`). |
| `--num_random_sets` | int | `1` | Generate N numbered random output files (requires `--method random` when >1). |

## Random Pair Set Generation

Use `--method random --num_random_sets N` to generate multiple random species-group files in one command.

Example:

```bash
esl-psc pairs \
  --tree_file /path/tree.nwk \
  --species_pheno_path /path/phenotypes.csv \
  --output_path /path/random_pairs/pairs.txt \
  --method random \
  --num_random_sets 100
```

This writes:
- `/path/random_pairs/pairs_001.txt`
- `/path/random_pairs/pairs_002.txt`
- ...
- `/path/random_pairs/pairs_100.txt`

If `--output_path` is given as a directory path ending in `/`, default names are used:
- `auto_pairs_groups_001.txt`, `auto_pairs_groups_002.txt`, ...

Seed behavior:
- If no `--seed` is provided: each set is random and non-reproducible.
- If `--seed S` is provided with multiple sets: runs use seeds `S, S+1, S+2, ...`.

## Method Notes

- `simple_deterministic` (or alias `default`): use the initial deterministic adjacent-transition pairing behavior.
- `shortest`: prioritize shortest phylogenetic distance tie-breaks.
- `longest`: tie-break by longer sequence evidence (requires `--alignments_dir`).
- `contrast`: in continuous mode, prefers larger trait contrast.
- `composite`: combines trait contrast, distance, and sequence support (requires `--alignments_dir`).
- `random`: random tie-breaking; supports `--seed` and `--num_random_sets`.
- `pct_contrast`: local percent-contrast selector for continuous positive phenotypes only; does not support alternates.

## Examples

Continuous thresholds:

```bash
esl-psc pairs \
  --tree_file /path/tree.nwk \
  --species_pheno_path /path/phenotypes.csv \
  --output_path /path/species_groups.txt \
  --method contrast \
  --lower_threshold 2.0 \
  --upper_threshold 4.0
```

Quantile-tail thresholding:

```bash
esl-psc pairs \
  --tree_file /path/tree.nwk \
  --species_pheno_path /path/phenotypes.csv \
  --output_path /path/species_groups.txt \
  --method contrast \
  --quantile_tails_pct 10
```

Percent-contrast local selector:

```bash
esl-psc pairs \
  --tree_file /path/tree.nwk \
  --species_pheno_path /path/phenotypes.csv \
  --output_path /path/species_groups.txt \
  --method pct_contrast \
  --min_pct_diff 25
```

## Validation and Common Errors

- `--num_random_sets > 1` requires `--method random`.
- `--method longest` and `--method composite` require `--alignments_dir`.
- `--method pct_contrast` requires continuous phenotype values that are strictly positive.
- At least 2 valid contrast pairs are required to produce output.
