# Evolutionary Sparse Learning with Paired Species Contrast (ESL-PSC) #

[![Release](https://img.shields.io/github/v/release/John-Allard/ESL-PSC?label=release)](https://github.com/John-Allard/ESL-PSC/releases/latest)
[![GUI builds](https://img.shields.io/github/actions/workflow/status/John-Allard/ESL-PSC/nuitka-build.yml?label=GUI%20builds)](https://github.com/John-Allard/ESL-PSC/actions/workflows/nuitka-build.yml)
[![Package artifacts](https://img.shields.io/github/actions/workflow/status/John-Allard/ESL-PSC/package-managers.yml?label=package%20artifacts)](https://github.com/John-Allard/ESL-PSC/actions/workflows/package-managers.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Table of Contents ##

1. [Description](#description)
2. [Quick Start](#quick-start)
3. [Installation and Dependencies](#installation-and-dependencies)
4. [Graphical User Interface](#graphical-user-interface)
5. [Command Line Usage](#command-line-usage)
   - [Auto Pair Selection CLI](#auto-pair-selection-cli)
   - [Site Counter CLI](#site-counter-cli)
6. [Additional CLI Tools](#additional-cli-tools)
7. [Input Data](#input-data)
8. [Output Data](#output-data)
9. [Additional Options and Parameters](#additional-options-and-parameters)
10. [Included Data](#included-data)
11. [Demo](#demo)
12. [Troubleshooting](#troubleshooting)
13. [Citation](#citation)

## Description ##
This repository provides tools for analyzing signatures of molecular convergence in a multiple sequence alignment using Evolutionary Sparse Learning with Paired Species Contrast (ESL-PSC). The main CLI, `esl-psc`, takes alignment files and a defined set of species contrast pairs as inputs, preprocesses input data and builds sparse group lasso models to explain the trait of interest. See [Allard et al., 2025, *Nature Communications*](https://www.nature.com/articles/s41467-025-58428-8) for methodological details. 

![flow chart](./images/ESL_PSC_flowchart_image.png)

## Quick Start ##

If you want to get started quickly, the GUI is the easiest way to run ESL-PSC without needing to assemble a CLI command manually. Start with the beginner guide here:

- [`docs/gui-quickstart.md`](docs/gui-quickstart.md)

It shows which GUI download to choose for Mac, Windows, or Ubuntu/Debian Linux, and walks through a first analysis with screenshots. If you prefer the terminal, the command-line tools are described below.

## Installation and Dependencies ##

### Install options

Prepackaged software is available for download.

#### 1. GUI app downloads (macOS/Windows/Linux)

Download the GUI app for your platform from [GitHub Releases](../../releases/latest):



#### 2. Toolkit archive (Linux/macOS/Windows)

The Toolkit provides the CLI for ESL-PSC. Download the toolkit archive for your platform from [GitHub Releases](../../releases/latest):



Then:

1. Extract archive.
2. Use Python 3.10+ to install toolkit Python dependencies:
   `python -m pip install -r requirements-toolkit.txt`
3. Run:
   - `./bin/esl-psc --help` (Linux/macOS)
   - `.\bin\esl-psc.exe --help` (Windows)

#### 3. Build/install from source (most reliable on Linux)

This avoids prebuilt binary issues on some Linux systems.

1. Install Rust toolchain.
2. Install CLI:
   `cargo install --path esl_psc_rs --root /your/install/prefix`
3. Install Python dependencies:
   `python -m pip install -r requirements-toolkit.txt`
4. Set Python linkage for utility subcommands:
   - `ESL_PSC_PYTHON` should point to your Python executable
   - `ESL_PSC_PYTHONPATH` should point to a directory containing `esl_psc_cli` and `gui` modules (the repo root works)

  

## Graphical User Interface ##

If you are new to ESL-PSC and want the simplest path, start with the guide in [`docs/gui-quickstart.md`](docs/gui-quickstart.md).


We introduced a GUI that wraps all of ESL-PSC's functionality in an easy-to-use wizard-style app:

Step-by-step workflow guiding you through:
  1. **Input selection** – Select files and directories for each required and optional input 
  2. **Analysis parameters** – adjust settings within allowed ranges using controls with helpful tooltips
  3. **Command preview** – view the exact CLI command that will be executed. This can be copied to the clipboard and run from the command line if you prefer.
  4. **Run** – Execute the analysis and view the terminal output within the GUI. When SPS plots are generated, they remain closed by default – use the **Show SPS Plot** button to open them. Running the same command directly from the CLI will display the plot automatically.
  5. **View results** – View gene rankings and examine convergent sites in an interactive alignment viewer.

Save your configuration in a file and re-load it later.

Now compatible with Windows, Mac, and Linux.

![ESL-PSC GUI](./images/ESL-PSC_GUI.png)

#### Interactive Tree Viewer & Automatic Contrast Pair Selection

The GUI now features an interactive phylogenetic tree viewer that lets you:

* Load a Newick tree and a species phenotype file to visualize convergent and non-convergent clades.
* Assign convergent vs. non-convergent phenotypes by right-clicking species names.
* Visualize phenotypes with colors:
  - Binary values (-1/1) color species labels red/blue.
  - Continuous float values are supported and colored by a red→blue gradient (low→high).
  - Species without phenotype data remain black.
  - Contrast pair coloring takes precedence over phenotype coloring to keep analysis pairs clear.
* Assign species to contrast pairs which are then marked on the tree view.
* Automatically generate sensible convergent/control contrast pairs with a single click.
* Save / load phenotype assignments and the species groups file to use in the analysis.  No need to create text or CSV files manually. 
* Export an SVG graphic of the annotated tree graphic with pairs selected to keep track of and display your experimental design.

Simply press **“Create a Species Groups File Using a Newick Tree”** on the first page, then use the toolbar at the top of the viewer.

![Tree Viewer](./images/tree_viewer.png)

#### Continuous Phenotype Support

We expanded GUI support for continuous (numeric) phenotypes across the Tree Viewer, analysis options, Site Viewer, and plots:

- Tree Viewer (continuous mode)
  - Load a phenotype file with float values. Species labels are colored by percentile using a Viridis gradient and the numeric value is shown next to each species name.
  - Automatic contrast pair selection supports two continuous workflows:
    - Thresholded binarization: prompt for lower/upper thresholds (with histogram + quantile tails) used only by the auto-selection algorithm.
    - Positive-trait percent contrast: require a minimum percent difference (`upper` vs `lower`) and greedily select closest non-overlapping pairs first.
  - For the positive-trait percent-contrast mode, the GUI can preview a threshold sweep plot (threshold vs selected pair count) before applying a threshold.
  - Tie-breaking options are available when multiple sibling choices are valid: Longest sequence, Shortest distance, Max trait contrast, Composite best, Random, or Simple deterministic.
  - Exported SVGs include the continuous colorbar and low/high labels. In continuous mode, manual binary tools such as Invert Phenotype and Set All to Non‑convergent are disabled by design.

- Analysis modes (Parameters page and CLI)
  - New GUI toggle: Use continuous phenotype values? This runs ESL-PSC with linear regression (ordinary least squares; sg_lasso_leastr) instead of binary logistic regression.
  - CLI flag: `--use_continuous_phenotypes` to enable continuous-response modeling. When active, binary SPS plots are disabled by default.
  - New plot option for continuous runs: a 2D density plot of true phenotype (X) vs SPS (Y). In the GUI this appears as Generate phenotype density plot when continuous mode is on. CLI flag: `--make_continuous_plot`.

- Site Viewer: phenotype-aware “Other Species” pane
  - Phenotype values (if provided) are displayed alongside species names, and you can sort “Other” species by Phenotype High → Low or Phenotype Low → High in addition to Alphabetical and MSA Order.

- Plots for predictions
  - Binary phenotypes: violin or KDE SPS plots remain available when a binary phenotype file is supplied.
  - Continuous phenotypes: select the new Phenotype vs SPS density plot (GUI checkbox or `--make_continuous_plot`).

#### Local Contrast Selector for Continuous Traits

- Tree Viewer (continuous mode)
  - Added a local percent-contrast selector for positive-valued continuous traits. This mode picks closest non-overlapping pairs that pass a user-chosen minimum percent-difference threshold.
  - Includes a threshold sweep preview (threshold vs selected pair count) to help choose a cutoff before applying pair selection.

- CLI auto pair selection
  - Added method `pct_contrast` with `--min_pct_diff` for the same local percent-contrast selection workflow.


### Running the GUI
1. Ensure the GUI dependencies are installed:
   ```bash
   pip install -r requirements-gui.txt  
   ```
   (You also need the core ESL-PSC requirements listed above.)
2. From the repository root folder, run:
   ```bash
   python -m gui.main
   ```
   The wizard window should appear. Fill in the pages and click **Run** on the final page to start the analysis.
   If you see an error about the `xcb` platform plugin needing
   `libxcb-cursor0`, install the system package first:
   ```bash
   sudo apt install libxcb-cursor0
   ```

### Additional dependencies
The GUI requires **PySide6** on top of other ESL-PSC requirements (`biopython`, `numpy`, `pandas`, `matplotlib`, `seaborn`). The `requirements-gui.txt` file lists everything you need.


Feedback on the GUI is welcome! Please open an issue on the [GitHub repository](https://github.com/John-Allard/ESL-PSC/issues) if you have any questions or suggestions.

### Stand-alone packaged applications ###

Pre-built GUI packages are available for macOS, Windows, and Linux on the [GitHub Releases page](../../releases/latest). We also publish a toolkit package for terminal use centered on `esl-psc` and its utility subcommands, including `pairs` and `site-counter`.

The toolkit package includes Rust binaries along with the remaining Python CLI modules and wrappers, and is intended to run with your system Python rather than bundling another Python runtime.
After extracting the toolkit, install dependencies with:

`python3 -m pip install -r requirements-toolkit.txt`

Toolkit release artifacts are versioned by platform as:
`esl-psc-toolkit-v<version>-<os>-<arch>.<tar.gz|zip>`, with companion
`.sha256` and `.manifest.json` files for integrity and metadata.

The conda recipe for packaging `esl-psc` is provided at
`packaging/conda/recipe/`.

Package-manager scaffolding is included for:

- conda-forge recipe rendering: `packaging/conda-forge/`
- Homebrew formula rendering: `packaging/homebrew/`
- apt/deb package building: `packaging/apt/`

#### macOS build (Apple-notarized)
1. Download the macOS GUI release asset.
2. Double-click the `.zip` to extract the `ESL-PSC.app` bundle.
3. Drag `ESL-PSC.app` to your **Applications** folder.
4. Open the app via Launchpad, Spotlight or Finder. Because the build is notarized, macOS should open it without additional warnings. If a confirmation dialog appears, click **Open**.

#### Windows build (currently unsigned)
1. Download the Windows GUI release asset.
2. Right-click the file and select **Extract All…** (or use your preferred unzip tool).
3. Inside the extracted folder, double-click `ESL-PSC.exe` to launch.
   • Windows SmartScreen will warn that the executable is unsigned. Click **More info** and then **Run anyway** to continue.

#### Linux build (Debian/Ubuntu)
1. Download the Linux GUI release asset `esl-psc-gui_<version>_amd64.deb`.
2. Install it:
   - GUI: double-click the `.deb` and click **Install**
   - Terminal: `sudo dpkg -i ./esl-psc-gui_<version>_amd64.deb`
3. Launch from your app menu (`ESL-PSC`) or run `esl-psc-gui`.

If you are on a non-Debian Linux distribution or prefer to run from source, use the installation options listed above.


## Command Line Usage ##
The main CLI, `esl-psc`, now uses a high-performance unified Rust implementation of the analysis pipeline. It preserves the core CLI behavior while substantially reducing runtime for many analyses.

For a concise overview of the unified command-line interface, configuration-file behavior, checkpointing, and utility subcommands, see:

- [`docs/commands/cli-reference.md`](docs/commands/cli-reference.md)

Quick links:

- Pair selection: [`docs/commands/pairs.md`](docs/commands/pairs.md)
- Site Counter: [`docs/commands/site-counter.md`](docs/commands/site-counter.md)

### Using a Configuration File with ESL-PSC ###

ESL-PSC scripts can utilize a configuration file to easily manage arguments that remain constant across multiple runs. The scripts will check for the presence of an esl_psc_config.txt file in the current working directory. If it exists, the function reads the arguments from the file and combines them with any additional command-line arguments provided when running the script. This allows you to keep common arguments in the configuration file, while providing run-specific arguments via the command line.

To use this feature:

1. An existing esl_psc_config.txt file is included in the directory with the ESL-PSC scripts. 

2. Add any desired arguments to the esl_psc_config.txt file, using the same format as you would when providing them on the command line. Each argument should be on a new line, followed by its corresponding value. For example: `--pheno_names C4 C3`

3. Save the `esl_psc_config.txt` file.

4. When running an ESL-PSC script, the scripts will automatically check for the presence of the esl_psc_config.txt file and incorporate its contents. Command-line arguments will override the values in the config file if both are provided for the same argument.

5. If the esl_psc_config.txt file is not found in the current working directory, the function will only parse command-line arguments.

## Additional CLI Tools ##

In addition to the main `esl-psc` analysis command, ESL-PSC includes two utility CLIs for common supporting tasks:

- `esl-psc pairs`: tools for building or assisting with species contrast pair selection workflows.
- `esl-psc site-counter`: tools for counting and summarizing convergent or phenotype-associated sites from alignments and analysis outputs.

See the command-specific documentation for details:

- [`docs/commands/pairs.md`](docs/commands/pairs.md)
- [`docs/commands/site-counter.md`](docs/commands/site-counter.md)

## Input Data ##

#### The main input files required for ESL-PSC are: ####

1. A directory of alignment files. These should be in **2-line FASTA format** and whose file names must have the file extension `.fas`. Each sequence should appear entirely on a single line below the line containing its identifier. If the sequence is split over multiple lines, the run will continue but a warning will be issued and the deletion canceler will convert the files to 2-line format, which may be slower. It is assumed that each seperate alignment file will be a different genomic component, such as a gene, a protein, an exon, a domain, etc. and each component will be treated as a "group" of sites in the analysis (see Methods in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8)). Use the argument `--alignments_dir` and give the full absolute path to the directory.

2. A species groups file.  This is a text file that contains a comma delimited list of species on each line. In the simplest case, one species identifier can be placed on each line. The first line must contain one or more species that possess the convergent trait under analysis, and the next line must contain one or more species that can serve as trait-negative controls for the species in the first line, such that the first two lines, and each subsequent pair of lines will define a contrast pair of species to use in the analysis (see [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8) for details on chosing contrast pairs for ESL-PSC analysis). When more than one species is given in a line, each of those species will be used in a seperate analysis, along with all combinations of other alternative speices.  Thus, the total number of species combinations can be calculated by the product of the number of species given on each line. In the analysis, species listed on the first line, and subsequent odd numbered lines, will be assigned a response value of 1, and the associated control species on the even numbered lines will be assigned a response value of -1. Use the argument `--species_groups_file` and give the full absolute path to the file.

#### Optional input files: ####

1. A species phenotype file.
   - CLI (analysis): a text file with each full species name followed by a comma and then a 1 or -1 for the true phenotype class to which that species belongs. A 1 typically refers to the convergent phenotype. If this file is not provided, the true phenotype will not be listed for each species prediction in the species_predictions output file. Use the argument `--species_pheno_path` and give the full absolute path to the file.
   - GUI Tree Viewer (visualization only): also accepts continuous float phenotype values. Species labels are colored on a red→blue gradient (low→high). Species without phenotype entries are shown in black. When contrast pairs are assigned, pair colors override phenotype colors for those species.

2. A directory of alignments to use for preditions. By default, any species in the input alignments that are not used in building any given model will be assigned a sequence prediction score (SPS) for that model, which will be included in the predictions output file. As an alternative, you can use a seperate directory of alignments for the predictions, however these still need to be fully aligned to any input species alignments or the predictions will be meaningless. Use the argument `--prediction_alignments_dir` and give the full absolute path to the directory.

3. Limited genes list. If you want to use a subset of the alignment files for model building without having to remove files from your alignments directory, provide a limited-genes list file. This is a plain-text file with one alignment file *name* per line – no directory paths. Each name must exactly match a `.fas` file in your alignments directory (including the `.fas` extension). Use the argument `--limited_genes_list` and give the full absolute path to this list file.

## Output Data ##

ESL-PSC generates two main types of output files: a Predictions File and a Gene Ranks File. These files will be placed in the directory specified by `--output_dir`.

#### Predictions File ####
The predictions file contains every prediction made by every model generated using every species combination in the analysis. Each line in the file lists the following information:

1. Species combination (an abrevaited list of the species used to train the model. for very large numbers of species, a name like combo_1 will be assigned instead for each combination)
2. Lambda1 (site sparsity hyperparameter)
3. Lambda2 (gene sparsity hyperparameter)
4. Penalty term (the constant term used to calculate the group penalty, see hyperparameters below for details)
5. Number of genes (the number of genes/proteins with sites included in the model)
6. Input Root Mean Squared Error (RMSE; this is referred to as the Model Fit Score (MFS) by [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8))
7. Species being predicted
8. Sequence Prediction Score (SPS) (a negative value indicates a prediction of the phenotype assigned a response value of -1 and a positive value indicates a prediction of opposite phenotype) 
9. True phenotype for the species (taken from the species_pheno_file if provided)

#### Gene Ranks File ####
The gene ranks file lists the genes (or proteins or other genomic components) used in the analysis, along with information about their rankings based on their model contributions. It is recommended to perform ontology enrichment tests and/or other follow-up analyses on the highest ranking ~1% of genetic elements. Each line in the file includes the following information:

1. Gene name (taken from the alignment file)
2. Number of species combinations in which the gene is ranked (i.e. number of combinations for which it recieved a non-zero GSS as part of any model)
3. Number of species combinations in which the gene is ranked among the top contributors (the percentage of genes to consider "top genes" by GSS in any model can be set using the `--top_rank_frac` argument.)
4. Highest ever Group Sparsity Score (GSS)
5. Best ever rank (the best ever rank, 1 being the best possible, recieved in any model)

## Additional Options and Parameters ##

Detailed command-line options and parameter descriptions are maintained in the CLI reference:

- [`docs/commands/cli-reference.md`](docs/commands/cli-reference.md#additional-options-and-parameters)

## Included Data ##

#### We have included two sample species_groups files for use in ESL-PSC alignments ####
1. `demo_data/photosynthesis/photo_single_LC_matrix_species_groups.txt` (the grass species with the closest contrast partners with the longest sequences (i.e. fewest gaps; used for photosynthesis analyses in [Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8)))
2. `demo_data/echolocation/orthomam_echo_species_groups.txt` (this can be used to reproduce the echolocation analyses using all 16 species combinations ([Allard et al., 2025](https://doi.org/10.1038/s41467-025-58428-8))) 

A species phenotype file for the grass species has also been included: `demo_data/photosynthesis/photo_species_phenotypes.txt`

#### We have included the protein sequence alignments used for ESL-PSC analyses by Allard et al. (2025). If you use these data, please cite these sources: ####


##### Grass chloroplast alignments which were used by Allard et al. (2025) were derived from:

Casola C, Li J. 2022. Beyond RuBisCO: convergent molecular evolution of multiple chloroplast genes in C4 plants. PeerJ 10:e12791 https://doi.org/10.7717/peerj.12791
More information regarding these alignments can be found in the supplemental information kindly provided online by these authors.

##### Mammalian protein sequence alignments for echolocators and their control species were derived from the OrthoMaM database:
https://orthomam.mbb.cnrs.fr/#

OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes Celine Scornavacca, Khalid Belkhir, Jimmy Lopez, Rémy Dernat, Frédéric Delsuc, Emmanuel J P Douzery, Vincent Ranwez Molecular Biology and Evolution, Volume 36, Issue 4, April 2019, Pages 861–862

## Troubleshooting ##

Problems with the inputs can cause segmentation fault errors in the ESL preprocess step. Here are some common causes of problems:
1. An incorrect file path. It is recommended to use absolute file paths. Dragging the file icon onto the terminal window is a good way to make sure the path is entered correctly.
2. Having an extra blank new line in one of the input files.
3. having a duplicate alignment file name.
4. It is very easy to miss adding a ".txt" or other extension to one of the files names in the run command.
5. Non-standard characters in the alignments, like "*" for stop codons will cause problems.
6. Sequences must be aligned to each other in each file.

## Demo ##

### Running the demo in the GUI ###
To run the C3/C4 demo using the graphical interface:
1. Launch the GUI (see the installation and launch instructions above).
2. On the Input page, choose:
   - `demo_data/photosynthesis/alignments`
   - `demo_data/photosynthesis/photo_single_LC_matrix_species_groups.txt`
   - `demo_data/photosynthesis/photo_species_phenotypes.txt`
3. Click **Next** through the wizard pages and press **Run Analysis** on the final page.
4. When the run finishes, press **Show SPS Plot** to view the violin plot, or **Show Top Gene Ranks** to view the most influential genes. You can double click a gene's name to open the protein sequence alignment in the convergent site viewer and examine the residues present in the input species and other species at the strongest convergent sites.

### Running the demo with the command line ###
You can run an ESL-PSC analysis of the C3/C4 trait with the included chloroplast data by following the steps below: 
1. Clone this repository
2. Make sure you have the dependencies installed (see [Installation and Dependncies](#installation-and-dependncies) above). You will need 
3. Navigate to the `ESL_PSC/` directory on your computer
4. Run this command from the ESL-PSC directory: `esl-psc --output_file_base_name demo_output --species_groups_file demo_data/photosynthesis/photo_single_LC_matrix_species_groups.txt --alignments_dir demo_data/photosynthesis/alignments/ --use_logspace --num_log_points 20 --cancel_only_partner --species_pheno_path demo_data/photosynthesis/photo_species_phenotypes.txt --make_sps_plot --pheno_names "C4" "C3"`
5. The expected run time is approximately 30 seconds on a standard desktop computer.
6. A set of violin plots depeicting the prediction scores for C3 and C4 species will be displayed on the screen. The gene ranks (`demo_output_gene_ranks.csv`) and species prediction (`demo_output_species_predictions.csv`) csv files will be found in the ESL_PSC directory
the plot should look like this:
![predictions violin plot](./images/demo_output_image.png)
7. See [Output Data](#output-data) above for descriptions of the fields in the output csv files.

### Reproducing the original manuscript results ###

If you need to reproduce the results from the initial ESL-PSC manuscript, use version `v0.1.0`.

## Citation ##
If you use this software in your research, please cite our paper:

John B. Allard, Sudip Sharma, Ravi Patel, Maxwell Sanderford, Koichiro Tamura, Slobodan Vucetic, Glenn S. Gerhard & Sudhir Kumar. Evolutionary sparse learning reveals the shared genetic basis of convergent traits. Nature Communications 16, 3217 (2025). https://doi.org/10.1038/s41467-025-58428-8
