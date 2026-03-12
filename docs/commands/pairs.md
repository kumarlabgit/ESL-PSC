# ESL-PSC Pairs CLI

`esl-psc pairs` generates an ESL-PSC species-groups file from a phylogenetic tree and a phenotype file.

## Quick Start

Show help:

`esl-psc pairs --help`

Minimum required arguments:

```bash
esl-psc pairs \
  --tree_file /path/to/tree.nwk \
  --species_pheno_path /path/to/phenotypes.txt \
  --output_path /path/to/species_groups.txt
```

## Common Options

- `--method` selection strategy:
  - `default`
  - `simple_deterministic`
  - `longest`
  - `shortest`
  - `contrast`
  - `composite`
  - `random`
  - `pct_contrast`
- `--num_alternates` number of alternate species per branch
- `--max_combinations` limit the total combinatorial expansion
- `--alignments_dir` optional alignment directory for sequence-aware strategies
- Continuous-threshold options:
  - `--lower_threshold`
  - `--upper_threshold`
  - `--quantile_tails_pct`
  - `--min_pct_diff`
- Randomization:
  - `--seed`
  - `--num_random_sets`

## Output

Writes a species-groups file that can be used directly with the main ESL-PSC pipeline:

`--species_groups_file /path/to/species_groups.txt`

## Full Reference

For the full option reference, see:

- `esl_psc_cli/auto_pairs_cli_README.md`
