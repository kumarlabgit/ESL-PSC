# ESL-PSC GUI Wizard

A graphical user interface for running ESL-PSC analyses with an intuitive wizard interface.

## Features Overview

### Core Workflow
- **Step-by-step wizard** for configuring ESL-PSC analyses
- **Intuitive file selection** with validation and feedback
- **Parameter configuration** with sensible defaults and helpful tooltips
- **Real-time progress tracking** during analysis
- **Command preview** - see and copy the exact CLI command before running
- **Results viewer** - examine gene rankings and convergent sites interactively
- **Configuration management** - save/load analysis configurations as JSON files
- **Cross-platform support** (Windows, macOS, Linux)

### Interactive Tree Viewer

Load and visualize phylogenetic trees with phenotype data:

- **Tree file support**: Newick (`.nwk`, `.newick`, `.tree`, `.tre`) and NEXUS (`.nexus`, `.nex`)
- **Phenotype visualization**:
  - Binary (-1/1) phenotypes color labels red/blue
  - Continuous float phenotypes render via Viridis gradient (low→high)
  - Pair-based coloring overrides phenotype colors for clarity
  - Species without phenotype data appear black
- **Interactive editing**:
  - Right-click species to assign convergent/non-convergent phenotypes
  - Manually assign species to contrast pairs
  - Drag to pan, scroll to zoom
- **Automatic contrast pair selection**:
  - One-click generation of sensible convergent/control pairs
  - **Continuous phenotype thresholds**: Interactive histogram with tail-percentage control
  - **Tie-breaking options**: Longest sequence, Shortest distance, Max trait contrast, Composite best, Random, or Default
  - **Alternate selection**: Generate multiple alternates per species for robust analysis
- **Export**: Save SVG graphics with annotations, colorbar, and pair labels
- **Tree validation**: Warns if basal split is not a bifurcation

### Continuous Phenotype Support

Full support for numeric (continuous) phenotypes:

- **Analysis mode**: Toggle "Use continuous phenotype values?" to run linear regression (OLS) instead of binary logistic
- **Tree viewer gradient**: Percentile-based Viridis coloring with numeric values displayed
- **Threshold picker**: Interactive histogram with quantile controls for binarizing continuous values during pair selection
- **Site viewer**: Display phenotype values alongside species names, sortable by phenotype magnitude
- **Continuous plots**: Phenotype vs SPS density plot (2D heatmap) instead of binary violin/KDE plots

### Site Counter

Quick convergence screening of alignments without running full ESL-PSC:

- **Rapid screening**: Scan thousands of genes in minutes
- **Convergence metrics**: Avg true convergence, control convergence, difference, CS sites
- **Multi-combo support**: Analyze multiple species combinations and rank by cross-combo consistency
- **Two-pair mode**: Generate all 2×2 pair combinations from species groups
- **Adjustable outgroup agreement**: Relax the requirement for 100% control-outgroup agreement (0–100%)
- **Outgroup options**:
  - **Single species outgroup** (traditional)
  - **Parsimony ancestral reconstruction** (NEW) - see [Ancestral Reconstruction](#ancestral-reconstruction-outgroup) below
- **Results viewer**: Sortable table with gene details, click to open site viewer
- **Export**: Save results as CSV for further analysis

### Site Viewer

Interactive alignment viewer for examining convergent sites:

- **Alignment display**: Color-coded residues with convergent sites highlighted
- **Species grouping**: Convergent, Control, and Other species in separate panes
- **Phenotype-aware sorting**: Sort "Other" species by phenotype (high→low, low→high) when available
- **Site navigation**: Jump to specific convergent sites
- **CCS highlighting**: Convergent sites marked with yellow background
- **Statistics panel**: Display convergent site metrics and sequence information

### Additional Features

- **Response matrix mode**: Use pre-computed response matrices instead of species groups
- **Output options**: Control which outputs are generated (gene ranks, predictions, plots, etc.)
- **Load existing output**: Browse and view results from previous runs
- **Checkpoint/resume**: Resume interrupted analyses (checkpoints stored per alignment set)
- **Multi-format alignment support**: `.fas`, `.fasta`, `.fa`, `.faa`

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

The wizard window will appear. Navigate through the pages and click **Run** on the final page to start the analysis.

## Ancestral Reconstruction Outgroup

### Overview

Site Counter now supports using a **parsimony-reconstructed ancestral sequence** as the outgroup instead of selecting a single species. This provides more robust convergence detection when no single species serves as an ideal outgroup, or when you want to infer the ancestral state at the MRCA of your analysis species.

### How It Works

1. **Tree Pruning**: For each alignment, the input tree is pruned to only include species present in that alignment
2. **MRCA Identification**: The MRCA of all analysis species is identified in both the full tree and each pruned tree
3. **Validation**: Alignments are automatically skipped if the MRCA is at the root of the pruned tree (no outgroup available)
4. **Fitch Parsimony**: The algorithm reconstructs ancestral sequences position-by-position using a two-pass approach:
   - **Downpass** (tips → root): Intersect child state sets if non-empty, else union
   - **Uppass** (root → tips): Assign states preferring parent state when valid
5. **Convergence Detection**: The reconstructed ancestral sequence is used as the outgroup for CCS detection

### Usage in GUI

1. **Start Site Counter** from the Input Page
2. In the **Outgroup Dialog**, select **"Use parsimony ancestral reconstruction"**
3. Click **Browse** to select a tree file (Newick or NEXUS format)
4. The tree is validated automatically:
   - ✓ **Green check**: Tree is valid for ancestral reconstruction
   - ✗ **Red X**: Tree is missing species or has other issues (see error message)
5. Continue with Site Counter as normal

### Tree File Requirements

Your tree must:
- Be in **Newick** or **NEXUS** format (`.nwk`, `.newick`, `.tree`, `.nexus`, etc.)
- Contain **all species** from your species groups file
- Have the MRCA of analysis species **not at the root** (outgroup species must exist in the tree)
- Use species names matching your alignment files

### Advantages

- **Robustness**: Not dependent on a single species with potentially unique substitutions
- **Accuracy**: Infers the actual ancestral state using all available outgroup information
- **Flexibility**: Handles alignments with different species subsets automatically
- **Quality Control**: Skips alignments without proper outgroup context (logged with reason)

### Example Workflow

**Scenario**: Studying photosynthesis loss across 4 plant species with a phylogeny including 2 outgroups.

1. Create species groups file with your 4 analysis species
2. Launch Site Counter, choose ancestral reconstruction
3. Select your tree file (e.g., `plants_tree.nwk`)
4. Site Counter will:
   - Prune tree to present species per alignment
   - Identify MRCA of analysis species
   - Skip alignments where MRCA is at root
   - Reconstruct ancestral sequence at MRCA
   - Use ancestral sequence for convergence detection
5. Results show genes with significant convergence, excluding genes lacking outgroup context

### Validation and Skipping

Alignments are **skipped** (not counted as errors) when:
- MRCA of analysis species is at the root of the pruned tree
- No outgroup species are present in that alignment
- Reconstruction fails for technical reasons

This ensures only alignments with appropriate outgroup context are analyzed.

### Technical Details

- **Algorithm**: Fitch parsimony (Fitch 1971, Systematic Zoology 20:406-416)
- **Implementation**: `gui/core/ancestral_reconstruction.py`
- **Tests**: Comprehensive test suite in `tests/test_ancestral_reconstruction.py`
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
