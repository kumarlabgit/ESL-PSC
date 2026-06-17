# ESL-PSC Site Counter CLI

`esl-psc site-counter` scans alignments and reports convergent-change site metrics.

## Quick Start

Show help:

`esl-psc site-counter --help`

Basic run:

```bash
esl-psc site-counter \
  --alignments_dir /path/to/alignments \
  --species_groups_file /path/to/species_groups.txt \
  --outgroup_species Species_Name \
  --output_path /path/to/site_counter_output.csv
```

## Tree-Based Ancestral Mode

Use `--tree_file` instead of a fixed outgroup to reconstruct ancestral states at the MRCA:

```bash
esl-psc site-counter \
  --alignments_dir /path/to/alignments \
  --species_groups_file /path/to/species_groups.txt \
  --tree_file /path/to/tree.nwk \
  --output_path /path/to/site_counter_output.csv
```

Related flags:

- `--require_unambiguous_mrca`
- `--compute_mrca_representative`
- `--two_pair_combos`
- `--min_out_ctrl_agreement`
- `--top_frac`

## Output

Writes site-level results to the CSV file given by `--output_path`.
