# ESL-PSC GUI Wizard

This directory contains the PySide6 desktop interface for ESL-PSC. The GUI
collects the same inputs used by the command-line program, writes an executable
command, runs the analysis, and opens the main output tables and alignment
viewer from within the application.

## Main windows and tools

### Analysis wizard

The main window is a `QWizard` with pages for input files, run parameters,
command review, execution, and output inspection. It supports:

- alignment directories using `.fas`, `.fasta`, `.fa`, or `.faa` files
- species groups files, response matrix inputs, and optional phenotype files
- separate prediction alignment directories
- JSON configuration save/load
- terminal output from the running `esl-psc` process
- access to gene ranks, species predictions, selected sites, and SPS plots

### Tree viewer

The tree viewer opens Newick (`.nwk`, `.newick`, `.tree`, `.tre`) and NEXUS
(`.nexus`, `.nex`) trees. It can load phenotype files, color species labels by
binary or continuous phenotype values, and write a species groups file for the
main ESL-PSC analysis.

Tree annotations include:

- binary phenotype labels colored by class
- continuous phenotype labels colored by percentile
- contrast-pair colors for species assigned to training pairs
- manually assigned phenotypes and contrast pairs
- exported SVG files with tree labels, pair labels, and continuous color bars

The automatic pair selector can choose contrast pairs from the tree. Available
tie-breaking modes include longest sequence, shortest distance, maximum trait
contrast, composite score, random choice, and the default deterministic method.
Continuous traits can be converted to high/low classes with explicit thresholds
or quantile-tail thresholds. Positive-valued continuous traits can also use the
local percent-contrast selector.

### Continuous phenotype options

Continuous phenotype values are used in three places:

- the tree viewer, where species labels are colored by phenotype value
- automatic pair selection, where thresholds or percent contrast rules define
  high and low training species
- ESL-PSC runs, where the Parameters page can request linear regression instead
  of binary logistic regression

For continuous runs, the GUI can request the phenotype-vs-SPS density plot
instead of binary SPS violin or KDE plots.

### Site Counter

The Site Counter dialog runs `esl-psc site-counter` from the GUI. It accepts the
same species groups design used by ESL-PSC and reports convergence-count
statistics for each alignment. The dialog supports:

- fixed single-species outgroups
- parsimony-reconstructed ancestral outgroups from a tree file
- multiple species combinations
- two-pair combinations generated from species groups
- configurable control/outgroup agreement thresholds
- CSV export

### Site viewer

The site viewer opens an alignment and displays residues in separate panes for
the convergent, control, and other species. It can highlight convergent sites,
show CCS-compatible sites, and sort the other-species pane by phenotype when
phenotype annotations are available.

### Existing-output viewer

The output viewer can open result files from previous ESL-PSC runs, including
gene ranks, selected-site summaries, species predictions, and SVG plots.

## Installation

1. Ensure you have Python 3.8 or later installed
2. Install the required dependencies:
   ```bash
   pip install -r requirements-gui.txt
   ```
   (This includes PySide6, biopython, numpy, pandas, matplotlib, seaborn)

## Running the GUI

From the project root directory, run:

```bash
python -m gui.main
```

The wizard window will appear. Fill in the pages and click **Run** on the final
page to start the analysis.

## Ancestral Reconstruction Outgroup

Site Counter can use a parsimony-reconstructed ancestral sequence as the
outgroup instead of a named extant species. In this mode, the outgroup state is
estimated at the MRCA of the species included in the convergence analysis.

### Algorithm

1. For each alignment, the input tree is pruned to species present in that
   alignment.
2. The MRCA of all analysis species is identified in the full tree and in the
   pruned tree.
3. Alignments are skipped if the MRCA is at the root of the pruned tree, because
   no outgroup species remain for that alignment.
4. Ancestral residues are reconstructed position-by-position with Fitch
   parsimony:
   - downpass (tips to root): intersect child state sets if non-empty, otherwise
     take their union
   - uppass (root to tips): assign states while preferring the parent state when
     it is valid
5. The reconstructed ancestral sequence is used as the outgroup sequence for CCS
   detection.

### GUI use

1. Start Site Counter from the Input page.
2. In the outgroup dialog, select **Use parsimony ancestral reconstruction**.
3. Click **Browse** and select a Newick or NEXUS tree file.
4. The tree is validated automatically:
   - green check: tree is valid for ancestral reconstruction
   - red X: tree is missing species or has another validation error
5. Continue with Site Counter.

### Tree File Requirements

The tree must:
- Be in **Newick** or **NEXUS** format (`.nwk`, `.newick`, `.tree`, `.nexus`, etc.)
- Contain **all species** from your species groups file
- Have the MRCA of analysis species **not at the root** (outgroup species must exist in the tree)
- Use species names matching your alignment files

### Validation and Skipping

Alignments are **skipped** (not counted as errors) when:
- MRCA of analysis species is at the root of the pruned tree
- No outgroup species are present in that alignment
- Reconstruction fails for technical reasons

Skipped alignments are logged with the reason they were skipped.

### Technical Details

- **Algorithm**: Fitch parsimony (Fitch 1971, Systematic Zoology 20:406-416)
- **Implementation**: `gui/core/ancestral_reconstruction.py`
- **Tests**: `tests/test_ancestral_reconstruction.py`
- **Performance**:
  - **Compiled backend**: Performs both ancestral reconstruction and CCS detection natively. End-to-end runtime stays within a few milliseconds per gene (≈6 s for 20 k genes on test hardware).
  - **Python fallback**: Available when Rust is missing. Single-species outgroups are modestly slower (~0.26 s for 67 genes). Parsimony ancestor mode is much heavier (~30 s for 67 genes) because Python must parse the tree and reconstruct every alignment; use only when Rust cannot run.
  - **Caching**: Tree parsing is cached internally to reduce repeated work when Python fallback is engaged.

## Development

### Project Structure

- `gui/main.py` - Application entry point
- `gui/ui/` - User interface components
  - `gui/ui/pages/` - Wizard pages (InputPage, ParametersPage, RunPage, ResultsPage)
  - `gui/ui/widgets/` - Reusable UI widgets (TreeViewer, dialogs, file selectors, etc.)
- `gui/core/` - Business logic and data handling
  - `gui/core/fast_scan.py` - Site Counter implementation
  - `gui/core/ancestral_reconstruction.py` - Parsimony-based ancestral reconstruction
  - `gui/core/worker.py` - Background task execution
  - `gui/core/config.py` - Configuration management
- `gui/resources/` - Icons, styles, and other resources

### Testing

Run GUI-specific tests:
```bash
python3 tests/test_ancestral_manual.py
```

For full test suite (if pytest is installed):
```bash
pytest tests/test_ancestral_reconstruction.py -v
```
