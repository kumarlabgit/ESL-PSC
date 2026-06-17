## Command Line Usage ##
The ESL-PSC command line interface is provided by the unified `esl-psc` executable. Most users will run it either from a downloaded Toolkit archive or from a source build.

If you installed the Toolkit archive, run the executable from the extracted archive unless you have added it to your `PATH`:

- Linux/macOS: `./bin/esl-psc --help`
- Windows: `.\bin\esl-psc.exe --help`

If you built from source with Cargo, run the installed executable if it is on your `PATH`, or run the release binary from the Rust build directory:

- `esl-psc --help`
- `./esl_psc_rs/target/release/esl-psc --help`

The main analysis pipeline can be run either directly (`esl-psc ...`) or explicitly as `esl-psc run ...`. Utility functionality is available as subcommands including `esl-psc pairs` and `esl-psc site-counter`.

You can provide input parameters and options directly on the command line or through an `esl_psc_config.txt` configuration file.

Here is an example of how to run the main analysis command when `esl-psc` is on your `PATH`:

`esl-psc --output_file_base_name output_file_name --species_groups_file /path/to/species_groups_file --alignments_dir /path/to/alignments/dir --use_logspace --cancel_only_partner`

To list CLI options, run `esl-psc --help`.

Strict line-search acceptance is now the default solver behavior and matches the original ESL-PSC paper-era compatibility mode. Use `--enable_ec` only if you explicitly want the newer epsilon-comparison line-search acceptance.

See the repository [Demo](../../README.md#demo) for an example of a run command you can try with an included data set.

## Basic Input and Output Flags

The examples below use underscore-style flags. The Rust executable also accepts hyphenated aliases for these long options, such as `--alignments-dir`.

Required for a standard species-groups run:

- `--alignments_dir`: directory containing alignment FASTA files (`.fas`, `.fasta`, `.fa`, or `.faa`). Use a full absolute path unless you are intentionally running from a directory where relative paths are unambiguous.
- `--species_groups_file`: text file defining paired species contrasts. Provide exactly one input design source: `--species_groups_file`, `--response_dir`, or `--response_matrix_path`.
- `--output_file_base_name`: base name used for output files, such as `<base>_gene_ranks.csv` and `<base>_species_predictions.csv`.

Common optional input flags:

- `--species_pheno_path`: species phenotype file with one species and phenotype value per line. Use this when you want true phenotypes in prediction outputs or prediction plots.
- `--prediction_alignments_dir`: separate alignment directory for prediction species. If omitted, predictions are made from species in the main alignment directory that were not used to train each model.
- `--limited_genes_list`: plain-text file with one alignment file name per line. Names must match files in `--alignments_dir`, including the FASTA extension.
- `--response_dir`: directory of precomputed response matrix files. Use this instead of `--species_groups_file` when you already have response matrices.
- `--response_matrix_path`: single response matrix file. Use this for a single-matrix run instead of a species-groups multimatrix run.

Common output and run-control flags:

- `--output_dir`: directory where output files and intermediate run folders are written. If omitted, ESL-PSC creates a timestamped output directory.
- `--no_pred_output`: skip the species predictions file. This is faster for proteome-scale analyses focused on candidate gene ranking.
- `--no_genes_output`: skip the gene ranks file.
- `--show_selected_sites`: add selected-site summaries to the gene ranks output and write `<base>_selected_sites.csv`.
- `--use_continuous_phenotypes`: use continuous phenotype values and linear regression rather than binary logistic regression.
- `--make_sps_plot`, `--make_sps_kde_plot`, `--make_continuous_plot`: generate prediction plots when phenotype information is available.

## Configuration Files

If a file named `esl_psc_config.txt` exists in the current working directory, `esl-psc` reads arguments from that file before adding arguments supplied directly on the command line. This is useful for settings that stay the same across many runs. Command-line arguments supplied at runtime override values from the configuration file when both are present.

Write arguments in the same form used on the command line, separated by whitespace or line breaks. For example:

```text
--pheno_names C4 C3
--use_logspace
--num_log_points 20
```

### Legacy Python Wrappers

Legacy Python wrappers are retained only for maintainability and historical reference. New analyses should use the unified `esl-psc` executable, which is substantially faster.

### Checkpointing & Resuming Interrupted Runs

ESL-PSC automatically checkpoints progress by default for runs with multiple species combinations. After each species combination finishes, it saves a compact record inside a `checkpoint/` folder within your `--output_dir`. When you rerun the *same* command, the script detects the checkpoint and resumes exactly where it left off—skipping finished combos, reusing existing gap-canceled alignments and preprocess directories, and continuing to checkpoint as it goes.

How it works in brief:

1. Command-line arguments from the current run are compared to those stored in `checkpoint/command.json`. Benign flags (e.g. `--make_sps_plot`) may differ; critical parameters must match.
2. If they match, ESL-PSC sets `--use_existing_alignments` and `--use_existing_preprocess` automatically and starts with the next unfinished combo.
3. If they do **not** match you will get a descriptive message and the run will stop, prompting you to either adjust your command or start fresh with `--force_from_beginning`.

Flags you can use:

* `--no_checkpoint` — Disable checkpointing entirely (In tests it takes minimal time so it is recommended to leave it on).
* `--force_from_beginning` — Delete any existing `checkpoint/` folder in `--output_dir` and start from scratch.

To resume after an unexpected interruption, simply rerun the original command.

## Options and Parameters ##

The following additional options and parameters can be specified when running ESL-PSC to fine-tune the analysis and control various aspects of the process. These options can be added as command line arguments or specified in the config file. 

Note that the word the word "gene" is used here to refer to the genomic components treated as "groups" of sites in the analysis, but these can of course be any aligned segments of biological sequences, e.g. proteins, exons, etc.

##### Hyperparameters:
* `--initial_lambda1`: Initial lambda 1 value (position sparsity parameter). Default = .01.
* `--final_lambda1`: Final lambda 1 value (position sparsity parameter). Default = .99.
* `--initial_lambda2`: Initial lambda 2 value (group sparsity parameter). Default = .01.
* `--final_lambda2`: Final lambda 2 value (group sparsity parameter). Default = .99.
* `--lambda_step`: The increment to increase the lambda values with each step. It is recommended to use a logspace (see options below) but in a linear gridsearch of sparsity hyperparameters, this controls the step between values.
* `--group_penalty_type`: Group penalty calculation type ("median", "sqrt", "linear", or "std"). Median will be used by default (see Methods in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8)). The setting called "std" (not recommended) will leave the standard lasso group penalties in place. The "linear" option will calculate group penalties as the number of variable sites in the gene's alignment plus a constant term that is the same across all genes. By default, this will be 1 for all genes, but it is also possible to use a range of different constant terms and repeat all model ensembles for each group penalty term. In order to do this, the initial, final and step can be set using the following three arguments.
* `--initial_gp_value`: Group penalty constant term initial value. See group_penalty_type above for explanation.
* `--final_gp_value`: Group penalty constant term final value. See group_penalty_type above for explanation.
* `--gp_step`: Group penalty constant term increment. The default is 6. See group_penalty_type above for explanation.
* `--num_log_points`: The number of values per sparsity hyperparameter (lambda1 and lambda2) in a logspace of values to test. Include the `--use_logspace` flag (see options below).
* `--pheno_names`: The names of the two phenotypes separated by a space, with the convergent phenotype coming first. by default "1" and "-1" will be used
* `--min_genes`: Minimum number of genes a model must have in order for that model to be included in the prediction scores plots. Default = 0.
* `--maxiter`: Maximum number of gradient-descent iterations performed by the optimizer for each model fit. Default = 100. 

##### Options:
* `--use_logspace`: *Recommended* Use a log space of points for lambda values instead of initial and final lambda values with a lambda step.
* `--use_existing_preprocess`: Attempt to reuse existing preprocess folders. If a
  required folder is missing it will be generated automatically and a warning
  will be printed.
* `--use_default_gp`: Don't replace group penalties (automatically set to True if the group_penalty_type is "std").
* `--output_dir`: Directory where all output will be stored. If not supplied, a folder named `<output_file_base_name>_<timestamp>` will be created one level above the ESL-PSC project directory. Intermediate folders like `preprocessed_data_and_models`, `gap-canceled_alignments` and `response_matrices` will be created inside this location as needed.
* `--keep_raw_output`: Don't delete the raw model output files for each run. The raw models can be found in the `preprocessed_data_and_models` directory within the directory specified by `--output_dir`. You can also set a new directory by using the `--esl_inputs_outputs_dir` argument, but note that any files ending in `.txt` will be cleared from this directory before each ESL-PSC run.
* `--show_selected_sites`: Output the top-scoring sites for each gene. When enabled, the gene ranks file gains a `num_selected_sites` column and a separate `<output_name>_selected_sites.csv` file lists each site with its PSS (Position Sparsity Score). **Positions in this file are 1-indexed for readability**, whereas positions in the raw model output remain 0-indexed.
* `--no_genes_output`: Don't output a gene ranks file. If only predictions output is desired, including the option will speed up the analysis.
* `--no_pred_output`: Don't output a species predictions file. If only gene ranks output is desired, including the option will significantly speed up the analysis.
* `--make_sps_plot`: Make a two-panel violin plot showing SPS density for each true phenotype. The left panel uses the lowest 5% of MFS models; the right panel shows species-averaged SPS over all models. SPS values beyond `-1` and `1` are clipped to those bounds in the violin display and labeled accordingly.
* `--make_sps_kde_plot`: Make a two-panel KDE plot showing SPS density for each true phenotype. The left panel uses the lowest 5% of MFS models; the right panel shows species-averaged SPS over all models (see Methods in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8))
* `--no_checkpoint`: Disable checkpointing entirely. No `checkpoint/` folder will be created and the run cannot be resumed later. Use this for quick exploratory runs where resumption is unnecessary.
* `--force_from_beginning`: When a checkpoint folder already exists in `--output_dir`, delete it and start the analysis from scratch. This is the recommended way to rerun an experiment after changing critical parameters or if the checkpoint has become corrupted.

##### Deletion Canceler Options:
* `--nix_full_deletions`: Don't create new files for fully canceled genes, i.e. if enough species are missing the entire alignment is excluded.
* `--cancel_only_partner`: Only cancel partner of any gap species at the site instead of eliminating the entire column.
* `--min_pairs`: The minimum number of pairs that must not have gaps or the whole site will be canceled.
* `--limited_genes_list`: Use only genes in this list. One file per line.
* `--canceled_alignments_dir`: Full path to the new alignments directory. Gap-canceled alignments for each species combo will be placed here. This may also be an existing folder of gap-canceled alignments for multimatrix ESL-PSC. 
* `--use_existing_alignments`: Use existing files in canceled_alignments_dir instead of running deletion canceling.

##### Multimatrix-specific Optional Arguments:
* `--top_rank_frac`: Fraction of genes to count as "top genes" for the purpose of rankings across multiple species combinations. A setting of 0 will result in counting the single highest ranked gene as a top gene. The default is 0.01 (1%).
* `--response_dir`: Folder with response matrices. You can provide this instead of using a species groups file. Any txt file in this folder is assumed to be a response matrix file.
* `--use_uncanceled_alignments`: Use the alignments_dir alignments for all matrices without doing gap canceling (not recommended).
* `--delete_preprocess`: Clear preprocess folders after each matrix run.
* `--make_null_models`: Make null response-flipped ESL-PSC models. Must have an even number of pairs. All balanced flippings of the response values will be generated for each combo and all will be run and aggregated to maximally decouple true convergences (see Methods in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8)). 
* `--make_pair_randomized_null_models`: Make null pair randomized ESL-PSC models. A copy of input deletion-canceled alignment will, for each variable site, be randomized such that the residues of each contrast pair will be either flipped or not and the ESL-PSC integration will be repeated for each one. The results are then aggregated for all (see Methods in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8)).
* `--num_randomized_alignments`: Number of pair-randomized alignments to make. Default is 10.
