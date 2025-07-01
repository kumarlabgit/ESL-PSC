"""Input-selection page of the ESL-PSC wizard."""
from PyQt6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QFrame, QRadioButton,
    QLabel, QButtonGroup, QFormLayout, QPushButton, QFileDialog, QMessageBox
)
import os

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
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
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
                "Directory containing alignment files in 2-line FASTA format. Each file must have the .fas extension. "
                "Each sequence must be entirely on a single line below its identifier. "
                "All sequences in a file must be aligned. Only standard amino acid and gap characters are allowed."
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
                "Subsequent lines follow the same pattern." 
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
        self.tree_btn = QPushButton("Help me pick species from a Newick tree")
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
        self.response_dir.path_changed.connect(
            lambda p: setattr(self.config, 'response_dir', p)
        )
        self.response_dir.setVisible(False)  # Start hidden
        self.input_files_layout.addWidget(self.response_dir)
        
        # Add input files frame to required layout
        req_layout.addWidget(self.input_files_frame)
        
        # Connect radio buttons to toggle visibility
        def update_input_visibility():
            use_response_dir = self.use_response_dir.isChecked()
            self.species_groups.setVisible(not use_response_dir)
            self.response_dir.setVisible(use_response_dir)
            
            # Update config to reflect the active input type
            if use_response_dir:
                # Clear species groups path when switching to response dir
                self.config.species_groups_file = ""
            else:
                # Clear response dir path when switching to species groups
                self.config.response_dir = ""
        
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
            description="Optional: comma-separated file with species phenotypes. "
            "First column is species ID, second column is phenotype value (1 or -1).\n"
            "This is *required* for species prediction analyses."
        )
        self.species_phenotypes.path_changed.connect(
            lambda p: setattr(self.config, 'species_phenotypes_file', p)
        )
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
            description="Optional: File containing list of alignment files. "
            "The analysis will be limited to only those alignments even if more are in your alignments directory."
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
        try:
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
                        if len(row) >= 2:
                            try:
                                phenos[row[0].strip()] = int(row[1])
                            except ValueError:
                                continue
            except Exception as exc:
                QMessageBox.warning(
                    self,
                    "Phenotypes Error",
                    f"Failed to parse phenotypes file:\n{exc}",
                )

        self._tree_window = TreeViewer(tree, phenotypes=phenos)
        self._tree_window.show()

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