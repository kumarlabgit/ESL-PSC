"""Input-selection page of the ESL-PSC wizard."""
from PySide6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QFrame, QRadioButton,
    QLabel, QButtonGroup, QFormLayout, QPushButton, QFileDialog, QMessageBox,
    QSizePolicy, QCheckBox
)
from PySide6.QtCore import Qt  # Needed for alignment

# New helper for viewing completed runs
from gui.ui.widgets.existing_output_viewer import select_and_show_existing_output
import os

from esl_psc_cli import esl_psc_functions as ecf

from gui.ui.widgets.file_selectors import FileSelector
from gui.ui.widgets.tree_viewer import TreeViewer
from Bio import Phylo
from .base_page import BaseWizardPage

class InputPage(BaseWizardPage):
    """Page for selecting input files and directories."""
    
    def __init__(self, config):
        """Initialize the input page."""
        super().__init__("Input Selection")
        self.config = config
        self.setSubTitle("Select the input files and directories for the analysis.")
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)

        # ─── Load Existing Output Button (top-left) ────────────────────
        load_prev_btn = QPushButton("Load and View Existing Output")
        load_prev_btn.setToolTip("Select an output folder from a completed ESL-PSC run to view its results.")
        load_prev_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        load_prev_btn.clicked.connect(lambda *_: select_and_show_existing_output(parent=self))
        container_layout.addWidget(load_prev_btn, alignment=Qt.AlignmentFlag.AlignRight)
        
        # Required inputs group
        req_group = QGroupBox("Required Inputs")
        req_layout = QVBoxLayout()
        
        # Input type selection
        input_type_group = QGroupBox("Input Type")
        input_type_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "Choose the type of input for your analysis. You can either provide a species groups file "
            "(recommended) or a directory containing pre-computed response matrices."
        )
        explanation.setWordWrap(True)
        input_type_layout.addWidget(explanation)
        
        # Radio buttons for input type
        self.input_type_group = QButtonGroup(self)
        
        self.use_species_groups = QRadioButton("Use species groups file (recommended)")
        self.use_species_groups.setChecked(True)
        self.use_species_groups.setToolTip(
            "Provide a text file with species groups, one per line. Each line should contain species "
            "names separated by tabs."
        )
        
        self.use_response_dir = QRadioButton("Use response matrix directory")
        self.use_response_dir.setToolTip(
            "Provide a directory containing pre-computed response matrix files. Use this only if you have "
            "already generated response matrices for your species combinations."
        )
        
        self.input_type_group.addButton(self.use_species_groups)
        self.input_type_group.addButton(self.use_response_dir)
        
        input_type_layout.addWidget(self.use_species_groups)
        input_type_layout.addWidget(self.use_response_dir)
        # The 'Use continuous phenotype values?' selector now lives on the Output Options page
        input_type_group.setLayout(input_type_layout)
        
        req_layout.addWidget(input_type_group)
        
        # Input files frame
        self.input_files_frame = QFrame()
        self.input_files_layout = QVBoxLayout(self.input_files_frame)
        self.input_files_layout.setContentsMargins(0, 10, 0, 0)  # Add some top margin
        
        # Alignment directory (always required)
        self.alignment_dir = FileSelector(
            "Alignment Directory:", 'directory',
            default_path=os.getcwd(),
            description=(
                "Directory containing alignment files in FASTA format. Each file must have the .fas extension. "
                "Sequences should be in 2-line FASTA format for best performance, but multi-line files are accepted and "
                "will be converted, which may be slower. All sequences in a file must be aligned and contain only standard "
                "amino acid and gap characters."
            )
        )
        self.alignment_dir.path_changed.connect(
            lambda p: setattr(self.config, 'alignments_dir', p)
        )
        self.input_files_layout.addWidget(self.alignment_dir)
        
        # Species groups file selector
        self.species_groups = FileSelector(
            "Species Groups File:", 'file',
            default_path=os.getcwd(),
            description=(
                "Text file that specifies contrast pairs of species to use in the analysis. "
                "Each pair goes on a pair of lines, with the first line containing a species with the convergent trait "
                "and optionally one or more close siblings with the trait. The second line contains one or more control species "
                "that do not have the trait but are close relatives of the species in the first line. "
                "Subsequent lines follow the same pattern. " 
                "Example:\n\n"
                "convergent_species1a, convergent_species1b\n"
                "control_species1a, control_species1b\n"
                "convergent_species2a, convergent_species2b\n"
                "control_species2\n"
                "..."
            )
        )
        self.species_groups.path_changed.connect(
            lambda p: setattr(self.config, 'species_groups_file', p)
        )
        self.input_files_layout.addWidget(self.species_groups)

        # Button to open a Newick tree viewer
        self.tree_btn = QPushButton("Create a Species Groups File Using a Newick Tree")
        self.tree_btn.setToolTip(
            "Load a Newick file in a tree viewer where you can select contrast pairs interactively, "
            "and save a species groups file to use in the analysis."
        )
        self.tree_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.tree_btn.clicked.connect(self.open_newick_tree)
        self.input_files_layout.addWidget(self.tree_btn)
        
        # Response directory selector (initially hidden)
        self.response_dir = FileSelector(
            "Response Matrix Directory:", 'directory',
            default_path=os.getcwd(),
            description=(
                "Directory containing pre-computed response matrix files. "
                "Only use this if you have already generated response matrices. "
                "Typically used for advanced analyses or when reusing previously computed matrices."
            )
        )
        self.response_dir.path_changed.connect(self._on_response_dir_changed)
        self.response_dir.setVisible(False)  # Start hidden
        self.input_files_layout.addWidget(self.response_dir)
        
        # Add input files frame to required layout
        req_layout.addWidget(self.input_files_frame)
        
        # Connect radio buttons to toggle visibility
        def update_input_visibility():
            use_response_dir = self.use_response_dir.isChecked()
            self.species_groups.setVisible(not use_response_dir)
            self.response_dir.setVisible(use_response_dir)
            # Show/hide the helper button for creating a species groups file
            self.tree_btn.setVisible(not use_response_dir)
            
            # Update config to reflect the active input type
            if use_response_dir:
                # Persist currently shown response dir in config
                if hasattr(self, 'response_dir'):
                    try:
                        self.config.response_dir = self.response_dir.get_path()
                    except Exception:
                        pass
                # Clear species groups path in config to avoid ambiguity for downstream code
                self.config.species_groups_file = ""
            else:
                # Persist currently shown species groups path in config
                if hasattr(self, 'species_groups'):
                    try:
                        self.config.species_groups_file = self.species_groups.get_path()
                    except Exception:
                        pass
                # Clear response dir path when switching to species groups
                self.config.response_dir = ""
                self.config.response_matrices_are_continuous = False
                self.config.use_continuous_phenotypes = False
            self._update_continuous_checkbox_visibility()
        
        self.use_species_groups.toggled.connect(update_input_visibility)
        self.use_response_dir.toggled.connect(update_input_visibility)
        
        req_group.setLayout(req_layout)
        container_layout.addWidget(req_group)
        
        # Optional inputs group
        opt_group = QGroupBox("Optional Inputs")
        opt_layout = QFormLayout()
        
        # Species phenotypes file
        self.species_phenotypes = FileSelector(
            "Species Phenotypes File:", 'file',
            default_path=os.getcwd(),
            description=(
                "Optional: comma-separated file with species phenotypes. "
                "First column is species ID, second column is phenotype value (-1/1 binary or float for continuous). "
                "If omitted, the predictions output will not include a true phenotype column."
            ),
        )
        self.species_phenotypes.path_changed.connect(self._on_pheno_path_changed)
        opt_layout.addRow(self.species_phenotypes)
        
        # Prediction alignments directory
        self.prediction_alignments = FileSelector(
            "Prediction Alignments Directory:", 'directory',
            default_path=os.getcwd(),
            description="Optional: Directory with alignments for prediction. If provided, will run prediction on these alignments."
        )
        self.prediction_alignments.path_changed.connect(
            lambda p: setattr(self.config, 'prediction_alignments_dir', p)
        )
        opt_layout.addRow(self.prediction_alignments)
        
        # Limited genes file
        self.limited_genes = FileSelector(
            "Limited Genes File:", 'file',
            default_path=os.getcwd(),
            description="Optional: Text file with one alignment file *name* per line (no directory paths). "
            "Each name must exactly match a `.fas` file in your alignments directory. The analysis will be limited to only those alignments even if more are present."
        )
        self.limited_genes.path_changed.connect(
            lambda p: setattr(self.config, 'limited_genes_file', p)
        )
        opt_layout.addRow(self.limited_genes)
        
        opt_group.setLayout(opt_layout)
        container_layout.addWidget(opt_group)

        
        
        # Add stretch to push everything to the top
        container_layout.addStretch()
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)

    def open_newick_tree(self):
        """Open a Newick file and display it in a tree viewer."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Newick File",
            os.getcwd(),
            "Newick Files (*.nwk *.newick *.tree *.txt);;All Files (*)",
        )
        if not path:
            return
        # Basic validation: ensure equal number of opening and closing parentheses
        try:
            with open(path, "r", errors="ignore") as _nf:
                newick_text = _nf.read()
            open_paren = newick_text.count("(")
            close_paren = newick_text.count(")")
            # Must have at least one pair of parentheses and counts must match
            if open_paren == 0 or close_paren == 0 or open_paren != close_paren:
                preview = newick_text[:100]
                raise ValueError(
                    "The file does not appear to be valid Newick: mismatched parentheses ("
                    f"{open_paren} '(' vs {close_paren} ')').\n\n"
                    f"File: {os.path.basename(path)}\nFirst 100 characters:\n{preview}"
                )

            tree = Phylo.read(path, "newick")
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to parse Newick file:\n{exc}",
            )
            return

        phenos = {}
        pheno_path = getattr(self.config, "species_phenotypes_file", "")
        if pheno_path and os.path.exists(pheno_path):
            try:
                import csv

                with open(pheno_path, newline="") as f:
                    reader = csv.reader(f)
                    for row in reader:
                        # Skip empty or too-short rows
                        if not row or len(row) < 2:
                            continue
                        name = row[0].strip()
                        val_str = row[1].strip()
                        if not name or not val_str:
                            continue
                        # Try to parse as float; ignore unparseable rows (e.g., header)
                        try:
                            val = float(val_str)
                        except ValueError:
                            continue
                        phenos[name] = val
            except Exception as exc:
                QMessageBox.warning(
                    self,
                    "Phenotypes Error",
                    f"Failed to parse phenotypes file (supports -1/1 or continuous floats):\n{exc}",
                )

        self._tree_window = TreeViewer(
            tree,
            phenotypes=phenos,
            on_pheno_changed=self._update_phenotype_file,
            on_groups_saved=self._update_groups_file,
            on_alignments_changed=self._update_alignment_dir,
            alignments_dir=getattr(self.config, 'alignments_dir', ''),
        )
        self._tree_window.show()

    # ──────────────────────────────────────────────────────────────────────────
    # Internal helpers
    # ──────────────────────────────────────────────────────────────────────────
    def _on_pheno_path_changed(self, path: str) -> None:
        """Set config path and detect phenotype type (binary vs continuous)."""
        setattr(self.config, 'species_phenotypes_file', path)
        # Default flags when no file
        self.config.species_pheno_is_binary = False
        self.config.species_pheno_is_continuous = False
        if not path or not os.path.exists(path):
            # Update visibility in case file was removed
            self._update_continuous_checkbox_visibility()
            # Inform Parameters page to refresh dependent UI immediately
            try:
                wiz = self.wizard()
                if wiz and hasattr(wiz, 'params_page') and hasattr(wiz.params_page, 'update_output_options_state'):
                    wiz.params_page.update_output_options_state()
            except Exception:
                pass
            return
        try:
            import csv
            has_value = False
            is_binary = True
            with open(path, newline="", errors="ignore") as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row or len(row) < 2:
                        continue
                    val_str = row[1].strip()
                    if not val_str:
                        continue
                    try:
                        val = float(val_str)
                    except ValueError:
                        # ignore header or malformed rows
                        continue
                    has_value = True
                    if val not in (-1.0, 1.0):
                        is_binary = False
                        break
            if has_value and is_binary:
                self.config.species_pheno_is_binary = True
                self.config.species_pheno_is_continuous = False
            elif has_value:
                self.config.species_pheno_is_binary = False
                self.config.species_pheno_is_continuous = True
        except Exception:
            # On error, leave flags False so downstream treats as absent for plots
            pass
        finally:
            self._update_continuous_checkbox_visibility()
            # Inform Parameters page to refresh dependent UI immediately
            try:
                wiz = self.wizard()
                if wiz and hasattr(wiz, 'params_page') and hasattr(wiz.params_page, 'update_output_options_state'):
                    wiz.params_page.update_output_options_state()
            except Exception:
                pass

    def _on_response_dir_changed(self, path: str) -> None:
        """Set response directory and detect continuous response matrices."""
        setattr(self.config, 'response_dir', path)
        self.config.response_matrices_are_continuous = False
        if not path or not os.path.isdir(path):
            self.config.use_continuous_phenotypes = False
            self._update_continuous_checkbox_visibility()
            # Inform Parameters page to refresh dependent UI immediately
            try:
                wiz = self.wizard()
                if wiz and hasattr(wiz, 'params_page') and hasattr(wiz.params_page, 'update_output_options_state'):
                    wiz.params_page.update_output_options_state()
            except Exception:
                pass
            return
        try:
            for fname in os.listdir(path):
                if not fname.endswith('.txt'):
                    continue
                full = os.path.join(path, fname)
                if ecf.response_matrix_is_continuous(full):
                    self.config.response_matrices_are_continuous = True
                    self.config.use_continuous_phenotypes = True
                    self.config.species_pheno_is_continuous = True
                    break
            else:
                self.config.use_continuous_phenotypes = False
                self.config.species_pheno_is_continuous = False
        except Exception:
            # Ignore detection errors; treat as binary
            self.config.use_continuous_phenotypes = False
            self.config.species_pheno_is_continuous = False
        finally:
            self._update_continuous_checkbox_visibility()
            # Inform Parameters page to refresh dependent UI immediately
            try:
                wiz = self.wizard()
                if wiz and hasattr(wiz, 'params_page') and hasattr(wiz.params_page, 'update_output_options_state'):
                    wiz.params_page.update_output_options_state()
            except Exception:
                pass

    def _update_continuous_checkbox_visibility(self) -> None:
        """Sync configuration flags for continuous mode; UI control is on Output Options page."""
        cont_detected = bool(
            getattr(self.config, 'species_pheno_is_continuous', False) or
            getattr(self.config, 'response_matrices_are_continuous', False)
        )
        if not cont_detected:
            # When no continuous inputs are present, ensure the setting is off
            self.config.use_continuous_phenotypes = False

    def _update_phenotype_file(self, path: str) -> None:
        """Update the phenotype file selector and config."""
        self.species_phenotypes.set_path(path)

    def _update_groups_file(self, path: str) -> None:
        self.species_groups.set_path(path)
        setattr(self.config, 'species_groups_file', path)

    def _update_alignment_dir(self, path: str) -> None:
        """Update the alignment directory selector and config."""
        self.alignment_dir.set_path(path)
        setattr(self.config, 'alignments_dir', path)

    # ──────────────────────────────────────────────────────────────────────────
    # Public helpers for wizard
    # ──────────────────────────────────────────────────────────────────────────
    def update_ui_from_config(self):
        """Synchronize widget states with values in self.config."""
        # Set paths for file/directory selectors
        if hasattr(self.config, 'alignments_dir'):
            self.alignment_dir.set_path(self.config.alignments_dir)
        if hasattr(self.config, 'species_groups_file'):
            self.species_groups.set_path(self.config.species_groups_file)
        if hasattr(self.config, 'response_dir'):
            self.response_dir.set_path(self.config.response_dir)
        if hasattr(self.config, 'species_phenotypes_file'):
            self.species_phenotypes.set_path(self.config.species_phenotypes_file)
        if hasattr(self.config, 'prediction_alignments_dir'):
            self.prediction_alignments.set_path(self.config.prediction_alignments_dir)
        if hasattr(self.config, 'limited_genes_file'):
            self.limited_genes.set_path(self.config.limited_genes_file)

        # Determine which input type is active based on populated paths
        use_resp = bool(self.config.response_dir)
        self.use_response_dir.setChecked(use_resp)
        self.use_species_groups.setChecked(not use_resp)
        # Sync the continuous toggle visibility and state
        self._update_continuous_checkbox_visibility()
