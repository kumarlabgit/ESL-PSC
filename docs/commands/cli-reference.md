## Command Line Usage ##
To use the ESL-PSC command line interface (CLI), run `esl-psc` with the necessary arguments and options. The main analysis pipeline can be run either directly (`esl-psc ...`) or explicitly as `esl-psc run ...`. Utility functionality is available as subcommands including `esl-psc pairs` and `esl-psc site-counter`.

You can provide the input parameters and options through the command line or by creating a configuration file called esl_psc_config.txt. When using a configuration file, provide one argument per line.

Here is an example of how to run the script:

`esl-psc --output_file_base_name output_file_name --species_groups_file /path/to/species_groups_file --alignments_dir /path/to/alignments/dir --use_logspace --cancel_only_partner`

To list CLI options, run `esl-psc --help`.

See [Demo](#demo) for an example of a run command you can try with an included data set.

Legacy Python wrappers are still available for maintainability and historical reference, but they are no longer kept in the repository root. They now live under `legacy/python_entrypoints/`.

If you need to run the legacy Python implementation directly, use module entry points such as:

- `python -m esl_psc_cli.esl_multimatrix`
- `python -m esl_psc_cli.auto_pairs_cli`
- `python -m esl_psc_cli.fast_scan_cli`
- `python -m esl_psc_cli.plot_cli`

### Auto Pair Selection CLI

Toolkit command: `esl-psc pairs --help`

Auto pair-selection docs:

- [`docs/commands/pairs.md`](docs/commands/pairs.md) (quick usage and examples)
- [`esl_psc_cli/auto_pairs_cli_README.md`](esl_psc_cli/auto_pairs_cli_README.md) (full option reference)

### Site Counter CLI

Toolkit command: `esl-psc site-counter --help`

Site Counter is integrated into the unified `esl-psc` binary and exposed through the `site-counter` subcommand. See [`esl_psc_rs/README.md`](esl_psc_rs/README.md) for unified CLI behavior and implementation details.

Site Counter docs:

- [`docs/commands/site-counter.md`](docs/commands/site-counter.md)

**Parsimony-based ancestral reconstruction** is available in Site Counter via the `--tree_file` option. Instead of specifying a single outgroup species, provide a phylogenetic tree and Site Counter will reconstruct the ancestral sequence at the MRCA (Most Recent Common Ancestor) of your analysis species for each alignment.

### NEW: Checkpointing & Resuming Interrupted Runs

ESL-PSC automatically checkpoints progress by default for runs with multiple species combinations. After each species combination finishes, it saves a compact record inside a `checkpoint/` folder within your `--output_dir`. When you rerun the *same* command, the script detects the checkpoint and resumes exactly where it left off—skipping finished combos, reusing existing gap-canceled alignments and preprocess directories, and continuing to checkpoint as it goes.

How it works in brief:

1. Command-line arguments from the current run are compared to those stored in `checkpoint/command.json`. Benign flags (e.g. `--make_sps_plot`) may differ; critical parameters must match.
2. If they match, ESL-PSC sets `--use_existing_alignments` and `--use_existing_preprocess` automatically and starts with the next unfinished combo.
3. If they do **not** match you will get a descriptive message and the run will stop, prompting you to either adjust your command or start fresh with `--force_from_beginning`.

Flags you can use:

* `--no_checkpoint` — Disable checkpointing entirely (In tests it takes minimal time so it is recommended to leave it on).
* `--force_from_beginning` — Delete any existing `checkpoint/` folder in `--output_dir` and start from scratch.

To resume after an unexpected interruption, simply rerun the original command.

