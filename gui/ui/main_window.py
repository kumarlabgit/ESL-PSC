"""
Main window for the ESL-PSC Wizard application.
"""
import os
from PyQt6.QtWidgets import (
    QMainWindow, QWizard, QWizardPage, QVBoxLayout, QLabel, QWidget,
    QGroupBox, QFormLayout, QScrollArea, QSizePolicy, QSpacerItem, QFrame, QStackedWidget,
    QDoubleSpinBox, QSpinBox, QComboBox, QCheckBox, QMessageBox, QButtonGroup,
    QTextEdit, QLineEdit, QPushButton, QProgressBar, QHBoxLayout, QApplication, QRadioButton,
    QFileDialog
)
from PyQt6.QtCore import Qt, pyqtSignal, QThreadPool
from PyQt6.QtGui import QFont, QTextCursor

from gui.core.worker import ESLWorker

from gui.core.config import ESLConfig
from gui.ui.widgets.file_selectors import FileSelector

class MainWindow(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        """Initialize the main window and set up the wizard."""
        super().__init__()
        print("MainWindow: Initializing...")
        
        try:
            # Set window properties
            self.setWindowTitle("ESL-PSC Wizard")
            self.setMinimumSize(1000, 750)  # Increased minimum size
            print("MainWindow: Window properties set")
            
            # Create a central widget
            print("MainWindow: Creating central widget...")
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            
            # Create main layout
            print("MainWindow: Creating main layout...")
            main_layout = QVBoxLayout(central_widget)
            main_layout.setContentsMargins(10, 10, 10, 10)
            main_layout.setSpacing(10)
            
            # Create the wizard
            print("MainWindow: Creating wizard...")
            self.wizard = ESLWizard()
            print("MainWindow: Wizard created")
            
            # Add wizard to layout
            main_layout.addWidget(self.wizard)
            
            # Set focus to the wizard
            self.wizard.setFocus()
            
            # Force update the layout
            self.update()
            self.repaint()
            
            print("MainWindow: Initialization complete")
            
        except Exception as e:
            print(f"MainWindow: Error during initialization: {str(e)}")
            import traceback
            print(traceback.format_exc())
            raise
            
    def showEvent(self, event):
        """Handle the show event to ensure proper window centering and visibility."""
        super().showEvent(event)
        
        # Center the window on screen
        screen = QApplication.primaryScreen().availableGeometry()
        size = self.size()
        if size.width() > screen.width() or size.height() > screen.height():
            self.resize(screen.size() * 0.8)  # 80% of screen size if too large
            
        # Center the window
        frame_geometry = self.frameGeometry()
        center_point = screen.center()
        frame_geometry.moveCenter(center_point)
        self.move(frame_geometry.topLeft())
        
        # Ensure the window is raised and activated
        self.raise_()
        self.activateWindow()
        print("MainWindow: Window shown and activated")


class ESLWizard(QWizard):
    """Multi-page wizard for configuring and running ESL-PSC analyses."""
    
    def __init__(self, parent=None):
        """Initialize the wizard."""
        print("ESLWizard: Initializing...")
        super().__init__(parent)
        print("ESLWizard: Parent initialized")
        
        try:
            # Set window properties first
            self.setObjectName("ESLWizard")
            self.setWindowTitle("ESL-PSC Analysis Wizard")
            self.setWizardStyle(QWizard.WizardStyle.ModernStyle)
            self.setOption(QWizard.WizardOption.HaveHelpButton, False)
            self.setOption(QWizard.WizardOption.HaveNextButtonOnLastPage, False)
            self.setOption(QWizard.WizardOption.HaveFinishButtonOnEarlyPages, False)  # No grayed-out finish button
            self.setOption(QWizard.WizardOption.NoBackButtonOnLastPage, True)
            self.setOption(QWizard.WizardOption.NoCancelButton, False)
            self.setMinimumSize(1000, 700)  # Slightly larger minimum size
            print("ESLWizard: Window properties set")
            
            # Apply stylesheet
            self.setStyleSheet("""
                QWizard {
                    background-color: #f5f5f5;
                }
                QWizardPage {
                    background-color: white;
                    padding: 20px;
                }
                QWizardPage > QLabel {
                    font-size: 14px;
                    margin-bottom: 10px;
                }
                QWizardPage > QGroupBox {
                    font-weight: bold;
                    border: 1px solid #ccc;
                    border-radius: 5px;
                    margin-top: 2ex;
                    padding: 10px;
                }
                QWizardPage > QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 3px 0 3px;
                }
            """)
            
            # Initialize configuration
            print("ESLWizard: Initializing config...")
            self.config = ESLConfig()
            print("ESLWizard: Config initialized")
            
            # Store page references
            print("ESLWizard: Creating pages...")
            self.input_page = InputPage(self.config)
            self.params_page = ParametersPage(self.config)
            self.output_page = OutputPage(self.config)
            self.run_page = RunPage(self.config)
            
            # Add pages
            print("ESLWizard: Adding pages...")
            self.addPage(self.input_page)
            print("ESLWizard: Added InputPage")
            self.addPage(self.params_page)
            print("ESLWizard: Added ParametersPage")
            self.addPage(self.output_page)
            print("ESLWizard: Added OutputPage")
            self.addPage(self.run_page)
            print("ESLWizard: Added RunPage")
            
            # Connect signals
            self.currentIdChanged.connect(self.on_current_id_changed)
            
            print("ESLWizard: Initialization complete")
        except Exception as e:
            print(f"ESLWizard: Error during initialization: {str(e)}")
            import traceback
            print(traceback.format_exc())
            raise
        
        # Connect signals
        self.currentIdChanged.connect(self.on_current_id_changed)
    
    def on_current_id_changed(self, page_id):
        """Handle page changes in the wizard."""
        current_page = self.currentPage()
        if hasattr(current_page, 'on_enter'):
            current_page.on_enter()
    
    def validateCurrentPage(self):
        """Validate the current page before allowing the user to proceed."""
        current_page = self.currentPage()
        
        # Check required fields on input page
        if current_page == self.input_page:
            if not self.config.alignment_dir:
                QMessageBox.warning(self, "Missing Required Field", 
                                  "Please select an alignment directory.")
                return False
            if not self.config.species_groups_file:
                QMessageBox.warning(self, "Missing Required Field",
                                  "Please select a species groups file.")
                return False
        
        return super().validateCurrentPage()


class BaseWizardPage(QWizardPage):
    """Base class for wizard pages with common functionality."""
    
    def __init__(self, title, parent=None):
        """Initialize the base wizard page."""
        print(f"BaseWizardPage: Initializing page '{title}'")
        super().__init__(parent)
        self.setTitle(title)
        print(f"BaseWizardPage: '{title}' - Setting up layout")
        # Only set the layout once
        if not self.layout():
            self.setLayout(QVBoxLayout())
        print(f"BaseWizardPage: '{title}' - Layout set")
    
    def on_enter(self):
        """Called when the page is entered. Can be overridden by subclasses."""
        pass


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
            "names separated by tabs. This is the recommended approach as it allows for more flexible analysis."
        )
        
        self.use_response_dir = QRadioButton("Use response matrix directory")
        self.use_response_dir.setToolTip(
            "Provide a directory containing pre-computed response matrix files. Use this only if you have "
            "already generated response matrices for your species groups."
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
                "Each file should represent a different genomic component (gene, protein, exon, domain, etc.)."
            )
        )
        self.alignment_dir.path_changed.connect(
            lambda p: setattr(self.config, 'alignment_dir', p)
        )
        self.input_files_layout.addWidget(self.alignment_dir)
        
        # Species groups file selector
        self.species_groups = FileSelector(
            "Species Groups File:", 'file',
            default_path=os.getcwd(),
            description=(
                "Text file with species groups, one per line. "
                "The first line contains species with the convergent trait (assigned value 1). "
                "The second line contains control species (assigned value -1). "
                "Each line should be a comma-separated list of species. Example:\n\n"
                "species1,species2,species3\n"
                "control1,control2,control3"
            )
        )
        self.species_groups.path_changed.connect(
            lambda p: setattr(self.config, 'species_groups_file', p)
        )
        self.input_files_layout.addWidget(self.species_groups)
        
        # Response directory selector (initially hidden)
        self.response_dir = FileSelector(
            "Response Matrix Directory:", 'directory',
            default_path=os.getcwd(),
            description=(
                "Directory containing pre-computed response matrix files. "
                "Only use this if you have already generated response matrices. "
                "Each file should be a text file with a .txt extension, where each line represents a response matrix. "
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
            description="Optional: Tab-delimited file with species phenotypes. First column is species ID, subsequent columns are phenotype values."
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
            description="Optional: File containing list of gene IDs to analyze. If not provided, all genes in the alignment directory will be used."
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


class ParametersPage(BaseWizardPage):
    """Page for setting analysis parameters."""
    
    def __init__(self, config, parent=None):
        super().__init__("Analysis Parameters", parent)
        self.config = config
        self.has_species_pheno = False  # Track if we have a species phenotype file
        
        # Store references to widgets that might be accessed after deletion
        self.widgets_initialized = False
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        
        # Create container widget for scroll area
        container = QWidget()
        self.container_layout = QVBoxLayout(container)
        self.container_layout.setContentsMargins(10, 10, 10, 10)
        self.container_layout.setSpacing(15)
        
        # Add container to scroll area
        scroll.setWidget(container)
        
        # Add scroll area to the page's layout
        layout = self.layout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(scroll)
        
        # ===== Hyperparameters Section =====
        hyper_group = QGroupBox("Hyperparameters")
        hyper_layout = QFormLayout()
        
        # Grid Type Selection
        grid_type_group = QHBoxLayout()
        self.grid_type = QButtonGroup()
        
        self.logspace_btn = QRadioButton("Logarithmic Grid (Recommended)")
        self.linear_btn = QRadioButton("Linear Grid")
        
        self.grid_type.addButton(self.logspace_btn, 0)
        self.grid_type.addButton(self.linear_btn, 1)
        
        self.logspace_btn.setChecked(True)
        self.use_logspace = True
        
        self.logspace_btn.toggled.connect(
            lambda checked: setattr(self.config, 'use_logspace', checked)
        )
        
        grid_type_group.addWidget(self.logspace_btn)
        grid_type_group.addWidget(self.linear_btn)
        grid_type_group.addStretch()
        
        hyper_layout.addRow("Grid Type:", grid_type_group)
        
        # Lambda 1 (Position Sparsity)
        lambda1_group = QHBoxLayout()
        
        self.initial_lambda1 = QDoubleSpinBox()
        self.initial_lambda1.setRange(0.0, 1.0)
        self.initial_lambda1.setSingleStep(0.01)
        self.initial_lambda1.setValue(0.01)
        self.initial_lambda1.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_lambda1', v)
        )
        
        self.final_lambda1 = QDoubleSpinBox()
        self.final_lambda1.setRange(0.0, 1.0)
        self.final_lambda1.setSingleStep(0.01)
        self.final_lambda1.setValue(0.99)
        self.final_lambda1.valueChanged.connect(
            lambda v: setattr(self.config, 'final_lambda1', v)
        )
        
        self.lambda1_step = QDoubleSpinBox()
        self.lambda1_step.setRange(0.01, 1.0)
        self.lambda1_step.setSingleStep(0.01)
        self.lambda1_step.setValue(0.1)
        self.lambda1_step.valueChanged.connect(
            lambda v: setattr(self.config, 'lambda1_step', v)
        )
        
        # Number of log points
        self.num_log_points = QSpinBox()
        self.num_log_points.setRange(5, 1000)
        self.num_log_points.setValue(20)
        self.num_log_points.valueChanged.connect(
            lambda v: setattr(self.config, 'num_log_points', v)
        )
        
        # Stacked widget to toggle between step and log points
        self.lambda1_step_widget = QWidget()
        step_layout = QHBoxLayout(self.lambda1_step_widget)
        step_layout.addWidget(QLabel("Step:"))
        step_layout.addWidget(self.lambda1_step)
        step_layout.setContentsMargins(0, 0, 0, 0)
        
        self.lambda1_log_widget = QWidget()
        log_layout = QHBoxLayout(self.lambda1_log_widget)
        log_layout.addWidget(QLabel("Points:"))
        log_layout.addWidget(self.num_log_points)
        log_layout.setContentsMargins(0, 0, 0, 0)
        
        self.lambda1_stack = QStackedWidget()
        self.lambda1_stack.addWidget(self.lambda1_step_widget)
        self.lambda1_stack.addWidget(self.lambda1_log_widget)
        self.lambda1_stack.setCurrentIndex(1 if self.logspace_btn.isChecked() else 0)
        
        # Connect grid type change to update the stack
        def update_grid_type():
            is_logspace = self.logspace_btn.isChecked()
            self.lambda1_stack.setCurrentIndex(1 if is_logspace else 0)
            self.lambda2_stack.setCurrentIndex(1 if is_logspace else 0)
            
        self.logspace_btn.toggled.connect(update_grid_type)
        self.linear_btn.toggled.connect(update_grid_type)
        
        lambda1_group.addWidget(QLabel("From:"))
        lambda1_group.addWidget(self.initial_lambda1)
        lambda1_group.addWidget(QLabel("To:"))
        lambda1_group.addWidget(self.final_lambda1)
        lambda1_group.addWidget(self.lambda1_stack)
        
        hyper_layout.addRow("Lambda 1 (Position Sparsity):", lambda1_group)
        
        # Lambda 2 (Group Sparsity)
        lambda2_group = QHBoxLayout()
        
        self.initial_lambda2 = QDoubleSpinBox()
        self.initial_lambda2.setRange(0.0, 1.0)
        self.initial_lambda2.setSingleStep(0.01)
        self.initial_lambda2.setValue(0.01)
        self.initial_lambda2.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_lambda2', v)
        )
        
        self.final_lambda2 = QDoubleSpinBox()
        self.final_lambda2.setRange(0.0, 1.0)
        self.final_lambda2.setSingleStep(0.01)
        self.final_lambda2.setValue(0.99)
        self.final_lambda2.valueChanged.connect(
            lambda v: setattr(self.config, 'final_lambda2', v)
        )
        
        self.lambda2_step = QDoubleSpinBox()
        self.lambda2_step.setRange(0.01, 1.0)
        self.lambda2_step.setSingleStep(0.01)
        self.lambda2_step.setValue(0.1)
        self.lambda2_step.valueChanged.connect(
            lambda v: setattr(self.config, 'lambda2_step', v)
        )
        
        # Stacked widget for Lambda 2 step/log points
        self.lambda2_step_widget = QWidget()
        step_layout2 = QHBoxLayout(self.lambda2_step_widget)
        step_layout2.addWidget(QLabel("Step:"))
        step_layout2.addWidget(self.lambda2_step)
        step_layout2.setContentsMargins(0, 0, 0, 0)
        
        self.lambda2_log_widget = QWidget()
        log_layout2 = QHBoxLayout(self.lambda2_log_widget)
        log_layout2.addWidget(QLabel("Points:"))
        log_layout2.addWidget(QLabel("20"))  # Just a label since it uses the same num_log_points
        log_layout2.setContentsMargins(0, 0, 0, 0)
        
        self.lambda2_stack = QStackedWidget()
        self.lambda2_stack.addWidget(self.lambda2_step_widget)
        self.lambda2_stack.addWidget(self.lambda2_log_widget)
        self.lambda2_stack.setCurrentIndex(1 if self.logspace_btn.isChecked() else 0)
        
        lambda2_group.addWidget(QLabel("From:"))
        lambda2_group.addWidget(self.initial_lambda2)
        lambda2_group.addWidget(QLabel("To:"))
        lambda2_group.addWidget(self.final_lambda2)
        lambda2_group.addWidget(self.lambda2_stack)
        
        hyper_layout.addRow("Lambda 2 (Group Sparsity):", lambda2_group)
        
        # Group Penalty Settings
        penalty_group = QVBoxLayout()
        
        # Penalty type selection
        type_group = QHBoxLayout()
        self.group_penalty_type = QComboBox()
        self.group_penalty_type.addItems(["median (Recommended)", "default", "sqrt", "linear"])
        self.group_penalty_type.setCurrentText("median (Recommended)")
        self.group_penalty_type.currentTextChanged.connect(
            lambda t: self._update_penalty_type(t.split(' ')[0])  # Remove " (Recommended)" from text
        )
        
        type_group.addWidget(QLabel("Type:"))
        type_group.addWidget(self.group_penalty_type)
        type_group.addStretch()
        
        penalty_group.addLayout(type_group)
        
        # Range settings (hidden by default for median)
        self.penalty_range_group = QWidget()
        range_layout = QHBoxLayout(self.penalty_range_group)
        
        self.initial_gp_value = QDoubleSpinBox()
        self.initial_gp_value.setRange(0.0, 100.0)
        self.initial_gp_value.setSingleStep(0.5)
        self.initial_gp_value.setValue(0.1)  # Default from CLI
        self.initial_gp_value.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_gp_value', v)
        )
        
        self.final_gp_value = QDoubleSpinBox()
        self.final_gp_value.setRange(0.0, 100.0)
        self.final_gp_value.setSingleStep(0.5)
        self.final_gp_value.setValue(1.0)  # Default from CLI
        self.final_gp_value.valueChanged.connect(
            lambda v: setattr(self.config, 'final_gp_value', v)
        )
        
        self.gp_step = QDoubleSpinBox()
        self.gp_step.setRange(0.1, 10.0)
        self.gp_step.setSingleStep(0.1)
        self.gp_step.setValue(0.1)  # Default from CLI
        self.gp_step.valueChanged.connect(
            lambda v: setattr(self.config, 'gp_step', v)
        )
        
        range_layout.addWidget(QLabel("From:"))
        range_layout.addWidget(self.initial_gp_value)
        range_layout.addWidget(QLabel("To:"))
        range_layout.addWidget(self.final_gp_value)
        range_layout.addWidget(QLabel("Step:"))
        range_layout.addWidget(self.gp_step)
        
        penalty_group.addWidget(self.penalty_range_group)
        
        # Initially hide the range group
        self.penalty_range_group.setVisible(False)
        
        hyper_layout.addRow("Group Penalty Settings:", penalty_group)
        
        # Phenotype Names
        pheno_group = QHBoxLayout()
        
        self.pheno_name1 = QLineEdit("Convergent")
        self.pheno_name1.setPlaceholderText("Positive phenotype name")
        self.pheno_name1.textChanged.connect(
            lambda t: setattr(self.config, 'pheno_name1', t)
        )
        
        self.pheno_name2 = QLineEdit("Control")
        self.pheno_name2.setPlaceholderText("Negative phenotype name")
        self.pheno_name2.textChanged.connect(
            lambda t: setattr(self.config, 'pheno_name2', t)
        )
        
        pheno_group.addWidget(QLabel("Positive:"))
        pheno_group.addWidget(self.pheno_name1)
        pheno_group.addWidget(QLabel("Negative:"))
        pheno_group.addWidget(self.pheno_name2)
        
        hyper_layout.addRow("Phenotype Names:", pheno_group)
        
        # Add Top rank fraction to Hyperparameters section
        self.top_rank_frac = QDoubleSpinBox()
        self.top_rank_frac.setRange(0.0001, 1.0)
        self.top_rank_frac.setSingleStep(0.01)
        self.top_rank_frac.setValue(0.01)
        self.top_rank_frac.setToolTip(
            "Fraction of top-ranked genes to consider for significance testing. Lower values make the test "
            "more stringent by focusing on the very top predictions."
        )
        self.top_rank_frac.valueChanged.connect(
            lambda v: setattr(self.config, 'top_rank_frac', v)
        )
        hyper_layout.addRow("Top rank fraction:", self.top_rank_frac)
        
        hyper_group.setLayout(hyper_layout)
        # Add hyper group to container
        self.container_layout.addWidget(hyper_group)
        
        # ===== Output Options =====
        output_group = QGroupBox("Output Options")
        output_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "Control which output files are generated. Including more outputs will increase "
            "computation time and disk usage."
        )
        explanation.setWordWrap(True)
        output_layout.addWidget(explanation)
        
        # Radio buttons for output options
        self.output_options_group = QButtonGroup(self)
        
        # Option 1: Gene ranks only (fastest)
        genes_only_layout = QHBoxLayout()
        self.genes_only_btn = QRadioButton("Gene ranks only (fastest)")
        self.genes_only_btn.setToolTip(
            "Only generate gene ranks output. This is the fastest option and is "
            "recommended for initial exploration."
        )
        self.genes_only_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_pred_output', checked)
        )
        self.output_options_group.addButton(self.genes_only_btn)
        genes_only_layout.addWidget(self.genes_only_btn)
        genes_only_layout.addStretch()
        output_layout.addLayout(genes_only_layout)
        
        # Option 2: Species predictions only (requires species phenotype file)
        preds_only_layout = QHBoxLayout()
        self.preds_only_btn = QRadioButton("Species predictions only (requires species phenotype file)")
        self.preds_only_btn.setToolTip(
            "Only generate species predictions output. This requires a species phenotype file "
            "to be specified on the input page."
        )
        self.preds_only_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_genes_output', checked)
        )
        self.output_options_group.addButton(self.preds_only_btn)
        preds_only_layout.addWidget(self.preds_only_btn)
        preds_only_layout.addStretch()
        output_layout.addLayout(preds_only_layout)
        
        # Option 3: Both outputs (slowest)
        both_outputs_layout = QHBoxLayout()
        self.both_outputs_btn = QRadioButton("Both outputs (slowest)")
        self.both_outputs_btn.setToolTip(
            "Generate both gene ranks and species predictions outputs. "
            "This is the most comprehensive but slowest option."
        )
        self.both_outputs_btn.setChecked(True)  # Default selection
        self.output_options_group.addButton(self.both_outputs_btn)
        both_outputs_layout.addWidget(self.both_outputs_btn)
        both_outputs_layout.addStretch()
        output_layout.addLayout(both_outputs_layout)
        
        # Output file base name
        output_file_layout = QHBoxLayout()
        output_file_layout.addWidget(QLabel("Output File Base Name:"))
        
        self.output_file_base_name = QLineEdit("esl_psc_results")
        self.output_file_base_name.setMaximumWidth(300)  # Make the field narrower
        self.output_file_base_name.textChanged.connect(
            lambda t: setattr(self.config, 'output_file_base_name', t)
        )
        # Set the default value in config
        self.config.output_file_base_name = "esl_psc_results"
        output_file_layout.addWidget(self.output_file_base_name, 1)
        output_layout.addLayout(output_file_layout)
        
        # Add vertical spacer for better separation
        output_layout.addSpacing(10)  # Add 10px spacing
        
        # Output directory
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(QLabel("Output Directory:"))
        
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setReadOnly(True)
        self.output_dir_edit.setPlaceholderText("Click Browse to select output directory")
        output_dir_layout.addWidget(self.output_dir_edit, 1)  # Allow expanding
        
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse_output_dir)
        output_dir_layout.addWidget(self.browse_btn)
        
        output_layout.addLayout(output_dir_layout)
        
        # Add vertical spacer for better separation
        output_layout.addSpacing(10)  # Add 10px spacing
        
        # Keep raw output
        self.keep_raw = QCheckBox("Keep raw output files")
        self.keep_raw.setToolTip("If checked, keep intermediate files generated during analysis.")
        self.keep_raw.stateChanged.connect(
            lambda s: setattr(self.config, 'keep_raw', s == 2)  # 2 is Qt.Checked
        )
        output_layout.addWidget(self.keep_raw)
        
        # Show selected sites
        self.show_selected_sites = QCheckBox("Show selected sites in output")
        self.show_selected_sites.setToolTip(
            "If checked, include a dictionary of all selected sites with their highest model score for every gene."
        )
        self.show_selected_sites.stateChanged.connect(
            lambda s: setattr(self.config, 'show_selected_sites', s == 2)
        )
        output_layout.addWidget(self.show_selected_sites)
        
        # Add vertical spacer for better separation
        output_layout.addSpacing(10)  # Add 10px spacing
        
        # SPS plot options (always visible but may be disabled)
        self.sps_plot_group = QGroupBox("Species Prediction Score (SPS) Plots")
        sps_plot_layout = QVBoxLayout()
        
        # Make SPS plot
        self.make_sps_plot = QCheckBox("Generate SPS density plots")
        self.make_sps_plot.setToolTip(
            "Create violin plots showing SPS density for each true phenotype."
        )
        self.make_sps_plot.stateChanged.connect(
            lambda s: setattr(self.config, 'make_sps_plot', s == 2)
        )
        sps_plot_layout.addWidget(self.make_sps_plot)
        
        # Make SPS KDE plot
        self.make_sps_kde_plot = QCheckBox("Generate SPS KDE plots")
        self.make_sps_kde_plot.setToolTip(
            "Create Kernel Density Estimate (KDE) plots showing SPS density for each true phenotype."
        )
        self.make_sps_kde_plot.stateChanged.connect(
            lambda s: setattr(self.config, 'make_sps_kde_plot', s == 2)
        )
        sps_plot_layout.addWidget(self.make_sps_kde_plot)
        
        self.sps_plot_group.setLayout(sps_plot_layout)
        output_layout.addWidget(self.sps_plot_group)
        
        # Connect output type changes to enable/disable SPS plot options
        def update_sps_plot_state():
            enable_sps = not self.genes_only_btn.isChecked()
            # Set the entire group box and its children enabled/disabled
            self.sps_plot_group.setEnabled(enable_sps)
            # Visually indicate the disabled state
            self.sps_plot_group.setStyleSheet(
                "QGroupBox:disabled { color: gray; }"
                "QCheckBox:disabled { color: gray; }"
            )
        
        self.genes_only_btn.toggled.connect(update_sps_plot_state)
        self.preds_only_btn.toggled.connect(update_sps_plot_state)
        self.both_outputs_btn.toggled.connect(update_sps_plot_state)
        update_sps_plot_state()  # Initial update
        
        output_group.setLayout(output_layout)
        # Add output group to container
        self.container_layout.addWidget(output_group)
        
        # ===== Deletion Canceler Options =====
        del_cancel_group = QGroupBox("Deletion Canceler Options")
        del_cancel_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns where gaps may be due to sequencing artifacts "
            "rather than true biological deletions. These settings control how the cancellation is performed."
        )
        explanation.setWordWrap(True)
        del_cancel_layout.addWidget(explanation)
        
        # Create form layout for the options
        form_layout = QFormLayout()
        
        # Nix full deletions
        self.nix_full_deletions = QCheckBox("Exclude fully deleted sites")
        self.nix_full_deletions.setToolTip(
            "If checked, sites that are fully deleted in any species will be excluded from analysis. "
            "This is equivalent to the --nix_full_deletions command line option."
        )
        self.nix_full_deletions.stateChanged.connect(
            lambda s: setattr(self.config, 'nix_full_deletions', s == 2)  # 2 is Qt.Checked
        )
        form_layout.addRow("Exclusion:", self.nix_full_deletions)
        
        # Cancel only partner
        self.cancel_only_partner = QCheckBox("Only cancel partner deletions")
        self.cancel_only_partner.setToolTip(
            "If checked, only cancel deletions that are part of a partner pair. "
            "If unchecked, all deletions will be canceled. This is equivalent to the --cancel_only_partner command line option."
        )
        self.cancel_only_partner.stateChanged.connect(
            lambda s: setattr(self.config, 'cancel_only_partner', s == 2)
        )
        form_layout.addRow("Deletion Mode:", self.cancel_only_partner)
        
        # Minimum aligned pairs
        min_pairs_layout = QHBoxLayout()
        self.min_pairs = QSpinBox()
        self.min_pairs.setRange(2, 100)  # Minimum is 2 as per requirements
        self.min_pairs.setValue(2)
        self.min_pairs.setToolTip(
            "Minimum number of aligned pairs required to consider a site. "
            "Sites with fewer aligned pairs will be excluded from analysis. "
            "This is equivalent to the --min_pairs command line option."
        )
        self.min_pairs.valueChanged.connect(
            lambda v: setattr(self.config, 'min_pairs', v)
        )
        min_pairs_layout.addWidget(self.min_pairs)
        min_pairs_layout.addStretch()
        form_layout.addRow("Minimum aligned pairs:", min_pairs_layout)
        
        # Add form layout to the main layout
        del_cancel_layout.addLayout(form_layout)
        
        del_cancel_group.setLayout(del_cancel_layout)
        self.container_layout.addWidget(del_cancel_group)
        
        # ===== Null Models Section =====
        null_models_group = QGroupBox("Null Models")
        null_models_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "Null models help assess the statistical significance of your results by comparing against "
            "randomized data. Select the type of null model to use:"
        )
        explanation.setWordWrap(True)
        null_models_layout.addWidget(explanation)
        
        # Radio buttons for null model selection
        self.null_models_group = QButtonGroup(self)
        
        # No null models (fastest)
        no_null_layout = QHBoxLayout()
        self.no_null_btn = QRadioButton("No null models")
        self.no_null_btn.setToolTip(
            "Do not generate any null models. This is the fastest option but provides "
            "no statistical significance estimates."
        )
        self.no_null_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_null_models', checked)
        )
        self.null_models_group.addButton(self.no_null_btn)
        no_null_layout.addWidget(self.no_null_btn)
        no_null_layout.addStretch()
        null_models_layout.addLayout(no_null_layout)
        
        # Response-flipped null models
        response_flip_layout = QHBoxLayout()
        self.response_flip_btn = QRadioButton("Response-flipped null models")
        self.response_flip_btn.setToolTip(
            "Generate null models by randomly flipping response values. "
            "This helps assess the significance of your results by comparing against random permutations "
            "of the phenotype assignments. Requires an even number of species pairs."
        )
        self.response_flip_btn.toggled.connect(
            lambda checked: setattr(self.config, 'make_null_models', checked)
        )
        self.null_models_group.addButton(self.response_flip_btn)
        response_flip_layout.addWidget(self.response_flip_btn)
        response_flip_layout.addStretch()
        null_models_layout.addLayout(response_flip_layout)
        
        # Pair-randomized null models
        pair_randomized_layout = QHBoxLayout()
        self.pair_randomized_btn = QRadioButton("Pair-randomized null models")
        self.pair_randomized_btn.setToolTip(
            "Generate null models by randomizing the alignment data while preserving the overall "
            "sequence composition. This controls for potential biases in the input data."
        )
        self.pair_randomized_btn.toggled.connect(
            lambda checked: setattr(self.config, 'make_pair_randomized_null_models', checked)
        )
        self.null_models_group.addButton(self.pair_randomized_btn)
        pair_randomized_layout.addWidget(self.pair_randomized_btn)
        pair_randomized_layout.addStretch()
        null_models_layout.addLayout(pair_randomized_layout)
        
        # Number of randomizations
        num_rand_layout = QHBoxLayout()
        self.num_rand = QSpinBox()
        self.num_rand.setRange(1, 1000)
        self.num_rand.setValue(10)
        self.num_rand.setToolTip(
            "Number of randomized alignments to generate for the null model. "
            "Higher values provide more accurate significance estimates but increase computation time. "
            "This is equivalent to the --num_randomized_alignments command line option."
        )
        self.num_rand.valueChanged.connect(
            lambda v: setattr(self.config, 'num_randomized_alignments', v)
        )
        num_rand_layout.addWidget(QLabel("Number of randomizations:"))
        num_rand_layout.addWidget(self.num_rand)
        num_rand_layout.addStretch()
        
        # Only show number of randomizations when a null model is selected
        self.num_rand_label = QLabel("Number of randomizations:")
        num_rand_layout = QHBoxLayout()
        num_rand_layout.addWidget(self.num_rand_label)
        num_rand_layout.addWidget(self.num_rand)
        num_rand_layout.addStretch()
        
        # Function to update visibility of num_rand based on selection
        def update_num_rand_visibility():
            show = not self.no_null_btn.isChecked()
            self.num_rand_label.setVisible(show)
            self.num_rand.setVisible(show)
        
        # Connect signals
        self.no_null_btn.toggled.connect(update_num_rand_visibility)
        self.response_flip_btn.toggled.connect(update_num_rand_visibility)
        self.pair_randomized_btn.toggled.connect(update_num_rand_visibility)
        
        # Initial update
        update_num_rand_visibility()
        
        null_models_layout.addLayout(num_rand_layout)
        
        # Set default selection
        self.no_null_btn.setChecked(True)
        
        null_models_group.setLayout(null_models_layout)
        self.container_layout.addWidget(null_models_group)
        
        # Mark widgets as initialized
        self.widgets_initialized = True
        
    def on_enter(self):
        """Called when the page is entered."""
        # Update our state based on the input page
        if hasattr(self.wizard(), 'input_page') and hasattr(self.wizard().input_page, 'species_phenotypes'):
            self.has_species_pheno = bool(self.wizard().input_page.species_phenotypes.get_path())
        
        # Update output options state when entering the page
        self.update_output_options_state()
        
        print("ParametersPage: Entering page")
        try:
            # Connect to update UI when input changes
            print("ParametersPage: Connecting to path_changed signal")
            self.wizard().input_page.species_phenotypes.path_changed.connect(
                self.update_output_options_state
            )
            # Initial update of output options state
            print("ParametersPage: Initial update of output options")
            self.update_output_options_state()
        except Exception as e:
            print(f"ParametersPage: Error in on_enter: {str(e)}")
            import traceback
            traceback.print_exc()
        
    def browse_output_dir(self):
        """Open a dialog to select the output directory."""
        try:
            # Get the current output directory or use the default
            current_dir = getattr(self.config, 'output_dir', os.getcwd())
            
            # Open directory dialog
            dir_path = QFileDialog.getExistingDirectory(
                self,
                "Select Output Directory",
                current_dir,
                QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks
            )
            
            if dir_path:  # User didn't cancel
                self.output_dir_edit.setText(dir_path)
                self.config.output_dir = dir_path
                
        except Exception as e:
            print(f"Error browsing for output directory: {e}")
            
    def update_output_options_state(self):
        """Update the state of output options based on whether a species phenotype file is provided."""
        if not hasattr(self, 'widgets_initialized') or not self.widgets_initialized:
            return
            
        try:
            # Check if widgets still exist
            if not hasattr(self, 'preds_only_btn') or not self.preds_only_btn:
                return
                
            # Update UI state based on our flag
            try:
                # Only update if state would change
                if self.preds_only_btn.isEnabled() != self.has_species_pheno:
                    self.preds_only_btn.setEnabled(self.has_species_pheno)
                    self.both_outputs_btn.setEnabled(self.has_species_pheno)
                
                # Update tooltips
                if self.has_species_pheno:
                    self.preds_only_btn.setToolTip(
                        "Only generate species predictions output. This requires a species phenotype file "
                        "to be specified on the input page."
                    )
                    self.both_outputs_btn.setToolTip(
                        "Generate both gene ranks and species predictions outputs. "
                        "This is the most comprehensive but slowest option."
                    )
                else:
                    self.preds_only_btn.setToolTip(
                        "<font color='red'>Requires species phenotype file to be specified on the input page.</font>"
                    )
                    self.both_outputs_btn.setToolTip(
                        "<font color='red'>Requires species phenotype file to be specified on the input page.</font>"
                    )
                
                # If species predictions were selected but no phenotype file, switch to gene ranks only
                if not self.has_species_pheno and (self.preds_only_btn.isChecked() or self.both_outputs_btn.isChecked()):
                    self.genes_only_btn.setChecked(True)
                    
            except RuntimeError:
                # Widgets have been deleted, ignore
                pass
                
        except Exception as e:
            print(f"Error updating output options state: {e}")
            # Ignore any other errors
            pass
            
    # ===== Deletion Canceler Options =====
    def _create_deletion_canceler_options(self):
        del_group = QGroupBox("Deletion Canceler Options")
        del_layout = QVBoxLayout()
        
        # Add the deletion canceler options to the container layout
        self.container_layout.addWidget(del_group)
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns where gaps may be due to sequencing artifacts "
            "rather than true biological deletions. These settings control how the cancellation is performed."
        )
        
        # Add the deletion canceler options to the container layout
        self.container_layout.addWidget(del_group)
        
        # Add stretch to push everything to the top
        self.container_layout.addStretch()
        
        # Set the container as the widget for the scroll area
        scroll.setWidget(container)
        
        # Add the scroll area to the page layout
        layout = QVBoxLayout(self)
        layout.addWidget(scroll)
        
        # Set the layout for the page
        self.setLayout(layout)
        
        # Mark widgets as initialized
        self.widgets_initialized = True
        
        # Connect signals last to avoid triggering updates during initialization
        if hasattr(self.wizard(), 'input_page'):
            if hasattr(self.wizard().input_page, 'species_phenotypes'):
                self.wizard().input_page.species_phenotypes.path_changed.connect(
                    self.update_output_options_state
                )
    
    def _create_deletion_canceler_options(self):
        """Create the deletion canceler options section."""
        del_group = QGroupBox("Deletion Canceler Options")
        del_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns where gaps may be due to sequencing artifacts "
            "rather than true biological deletions. These settings control how the cancellation is performed."
        )
        explanation.setWordWrap(True)
        del_layout.addWidget(explanation)
        
        # Add your deletion canceler options here
        
        del_group.setLayout(del_layout)
        self.container_layout.addWidget(del_group)
    
    def _update_penalty_type(self, penalty_type):
        """
        Update the UI based on the selected group penalty type.
        
        Args:
            penalty_type: The selected penalty type ('median', 'default', 'sqrt', or 'linear')
        """
        # Update the config
        self.config.group_penalty_type = penalty_type
        
        # Show/hide the range settings based on the penalty type
        if penalty_type == 'median':
            self.penalty_range_group.setVisible(False)
        else:
            self.penalty_range_group.setVisible(True)


class OutputPage(BaseWizardPage):
    """Page for reviewing configuration before running the analysis."""
    
    def __init__(self, config):
        """Initialize the output page."""
        super().__init__("Review Configuration")
        self.config = config
        self.setSubTitle("Review your configuration and command before running the analysis.")
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
        # Command preview section
        cmd_group = QGroupBox("ESL-PSC Command")
        cmd_layout = QVBoxLayout()
        
        # Command display with monospace font for better readability
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        self.cmd_display.setFont(QFont("Courier", 10))
        self.cmd_display.setPlaceholderText("The ESL-PSC command will be generated here...")
        self.cmd_display.setMinimumHeight(150)
        
        # Copy button
        self.copy_btn = QPushButton("Copy to Clipboard")
        self.copy_btn.clicked.connect(self.copy_command_to_clipboard)
        
        # Add widgets to command layout
        cmd_layout.addWidget(QLabel("This command will be executed:"))
        cmd_layout.addWidget(self.cmd_display)
        cmd_layout.addWidget(self.copy_btn, 0, Qt.AlignmentFlag.AlignRight)
        cmd_group.setLayout(cmd_layout)
        
        # Add command group to container
        container_layout.addWidget(cmd_group)
        
        # Output directory section
        output_dir_group = QGroupBox("Output Directory")
        output_dir_layout = QVBoxLayout()
        
        # Output directory display
        output_dir_widget = QWidget()
        output_dir_hbox = QHBoxLayout(output_dir_widget)
        output_dir_hbox.setContentsMargins(0, 0, 0, 0)
        
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setReadOnly(True)
        self.output_dir_edit.setStyleSheet("background-color: #f0f0f0; border: 1px solid #ccc; padding: 5px;")
        self.output_dir_edit.setToolTip("Output directory where results will be saved")
        
        # Set default output directory if not already set
        if not hasattr(self.config, 'output_dir') or not self.config.output_dir:
            self.config.output_dir = os.path.join(os.getcwd(), 'esl_psc_output')
        
        self.output_dir_edit.setText(self.config.output_dir)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.setMaximumWidth(80)
        browse_btn.clicked.connect(self.browse_output_dir)
        
        output_dir_hbox.addWidget(QLabel("Output Directory:"))
        output_dir_hbox.addWidget(self.output_dir_edit, 1)
        output_dir_hbox.addWidget(browse_btn)
        
        output_dir_layout.addWidget(output_dir_widget)
        output_dir_group.setLayout(output_dir_layout)
        
        # Configuration summary section
        summary_group = QGroupBox("Configuration Summary")
        summary_layout = QFormLayout()
        
        # Summary fields will be populated in on_enter()
        self.summary_labels = {}
        
        # Add some spacing between fields
        summary_layout.setVerticalSpacing(5)
        summary_group.setLayout(summary_layout)
        
        # Add output directory group to container before command group
        container_layout.insertWidget(0, output_dir_group)
        
        # Add summary group to container
        container_layout.addWidget(summary_group)
        
        # Add stretch to push everything to the top
        container_layout.addStretch()
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)
    
    def copy_command_to_clipboard(self):
        """Copy the current command to the system clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.cmd_display.toPlainText())
        
    def browse_output_dir(self):
        """Open a directory dialog to select the output directory."""
        current_dir = self.output_dir_edit.text() or os.getcwd()
        dir_path = QFileDialog.getExistingDirectory(
            self, 
            "Select Output Directory",
            current_dir,
            options=QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks
        )
        
        if dir_path:  # User didn't cancel
            self.config.output_dir = dir_path
            self.output_dir_edit.setText(dir_path)
            # Update the command display to reflect the new output directory
            self.on_enter()
        
    def on_enter(self):
        """Update the command preview and summary when the page is shown."""
        # Build the command
        cmd_parts = ["esl-psc"]
        
        # Required parameters
        if hasattr(self.config, 'alignment_dir') and self.config.alignment_dir:
            cmd_parts.append(f'--alignment-dir "{self.config.alignment_dir}"')
        if hasattr(self.config, 'species_groups_file') and self.config.species_groups_file:
            cmd_parts.append(f'--species-groups "{self.config.species_groups_file}"')
        
        # Optional parameters
        if hasattr(self.config, 'species_phenotypes_file') and self.config.species_phenotypes_file:
            cmd_parts.append(f'--species-phenotypes "{self.config.species_phenotypes_file}"')
        if hasattr(self.config, 'prediction_alignments_dir') and self.config.prediction_alignments_dir:
            cmd_parts.append(f'--prediction-alignments "{self.config.prediction_alignments_dir}"')
        if hasattr(self.config, 'limited_genes_file') and self.config.limited_genes_file:
            cmd_parts.append(f'--limited-genes "{self.config.limited_genes_file}"')
        if hasattr(self.config, 'output_dir') and self.config.output_dir:
            cmd_parts.append(f'--output-dir "{self.config.output_dir}"')
        
        # Hyperparameters
        if hasattr(self.config, 'initial_lambda1') and hasattr(self.config, 'final_lambda1') and hasattr(self.config, 'lambda1_step'):
            cmd_parts.append(f'--lambda1 {self.config.initial_lambda1}:{self.config.final_lambda1}:{self.config.lambda1_step}')
        
        if hasattr(self.config, 'initial_lambda2') and hasattr(self.config, 'final_lambda2') and hasattr(self.config, 'lambda2_step'):
            cmd_parts.append(f'--lambda2 {self.config.initial_lambda2}:{self.config.final_lambda2}:{self.config.lambda2_step}')
        
        # Group penalty settings
        if hasattr(self.config, 'group_penalty_type') and self.config.group_penalty_type != 'default':
            cmd_parts.append(f'--group-penalty {self.config.group_penalty_type}')
        
        if hasattr(self.config, 'initial_gp_value') and hasattr(self.config, 'final_gp_value') and hasattr(self.config, 'gp_step'):
            cmd_parts.append(f'--group-penalty-value {self.config.initial_gp_value}:{self.config.final_gp_value}:{self.config.gp_step}')
        
        # Logspace settings
        if hasattr(self.config, 'use_logspace') and self.config.use_logspace:
            cmd_parts.append(f'--logspace {self.config.num_log_points}')
        
        # Phenotype names
        if hasattr(self.config, 'pheno_name1') and hasattr(self.config, 'pheno_name2'):
            cmd_parts.append(f'--pheno-names {self.config.pheno_name1}:{self.config.pheno_name2}')
        
        # Minimum genes per model
        if hasattr(self.config, 'min_genes') and self.config.min_genes > 0:
            cmd_parts.append(f'--min-genes {self.config.min_genes}')
        
        # Deletion canceler options
        if hasattr(self.config, 'nix_full_deletions') and self.config.nix_full_deletions:
            cmd_parts.append('--nix-full-deletions')
        
        if hasattr(self.config, 'cancel_only_partner') and not self.config.cancel_only_partner:
            cmd_parts.append('--no-cancel-only-partner')
        
        if hasattr(self.config, 'min_pairs') and self.config.min_pairs > 1:
            cmd_parts.append(f'--min-pairs {self.config.min_pairs}')
        
        # Output options
        if hasattr(self.config, 'output_file_base_name') and self.config.output_file_base_name:
            cmd_parts.append(f'--output-file-base-name {self.config.output_file_base_name}')
        
        if hasattr(self.config, 'keep_raw_output') and self.config.keep_raw_output:
            cmd_parts.append('--keep-raw-output')
        
        if hasattr(self.config, 'show_selected_sites') and self.config.show_selected_sites:
            cmd_parts.append('--show-selected-sites')
        
        # Plot options
        if hasattr(self.config, 'make_sps_plot') and not self.config.make_sps_plot:
            cmd_parts.append('--no-sps-plot')
        
        if hasattr(self.config, 'make_sps_kde_plot') and self.config.make_sps_kde_plot:
            cmd_parts.append('--sps-kde-plot')
        
        # Multi-matrix options
        if hasattr(self.config, 'top_rank_frac') and self.config.top_rank_frac != 0.01:
            cmd_parts.append(f'--top-rank-frac {self.config.top_rank_frac}')
        
        if hasattr(self.config, 'response_dir') and self.config.response_dir:
            cmd_parts.append(f'--response-dir "{self.config.response_dir}"')
        
        # Null model options
        if hasattr(self.config, 'make_null_models') and self.config.make_null_models:
            cmd_parts.append('--make-null-models')
        
        if hasattr(self.config, 'make_pair_randomized_null_models') and self.config.make_pair_randomized_null_models:
            cmd_parts.append('--make-pair-randomized-null-models')
            if hasattr(self.config, 'num_randomized_alignments') and self.config.num_randomized_alignments != 10:
                cmd_parts.append(f'--num-randomized-alignments {self.config.num_randomized_alignments}')
        
        # Join command parts with spaces and proper line continuation
        cmd_str = r' \
  '.join(cmd_parts)
        
        # Display the command
        self.cmd_display.setPlainText(cmd_str)
        
        # Update the configuration summary
        self.update_summary()
    
    def update_summary(self):
        """Update the configuration summary section."""
        # Get the summary group box
        summary_group = self.findChild(QGroupBox, "Configuration Summary")
        if not summary_group:
            return
            
        # Get the layout and clear existing widgets
        layout = summary_group.layout()
        self.clear_layout(layout)
        
        # Add configuration items
        # Required parameters
        if hasattr(self.config, 'alignment_dir') and self.config.alignment_dir:
            self.add_summary_item(layout, "Alignment Directory:", self.config.alignment_dir)
        if hasattr(self.config, 'species_groups_file') and self.config.species_groups_file:
            self.add_summary_item(layout, "Species Groups File:", self.config.species_groups_file)
        
        # Optional input files
        if hasattr(self.config, 'species_phenotypes_file') and self.config.species_phenotypes_file:
            self.add_summary_item(layout, "Phenotypes File:", self.config.species_phenotypes_file)
        if hasattr(self.config, 'prediction_alignments_dir') and self.config.prediction_alignments_dir:
            self.add_summary_item(layout, "Prediction Alignments Dir:", self.config.prediction_alignments_dir)
        if hasattr(self.config, 'limited_genes_file') and self.config.limited_genes_file:
            self.add_summary_item(layout, "Limited Genes File:", self.config.limited_genes_file)
        
        # Hyperparameters
        if hasattr(self.config, 'initial_lambda1') and hasattr(self.config, 'final_lambda1') and hasattr(self.config, 'lambda1_step'):
            lambda1_str = f"{self.config.initial_lambda1} to {self.config.final_lambda1} (step: {self.config.lambda1_step})"
            self.add_summary_item(layout, "Lambda 1 Range:", lambda1_str)
        
        if hasattr(self.config, 'initial_lambda2') and hasattr(self.config, 'final_lambda2') and hasattr(self.config, 'lambda2_step'):
            lambda2_str = f"{self.config.initial_lambda2} to {self.config.final_lambda2} (step: {self.config.lambda2_step})"
            self.add_summary_item(layout, "Lambda 2 Range:", lambda2_str)
        
        # Group penalty settings
        if hasattr(self.config, 'group_penalty_type') and hasattr(self.config, 'initial_gp_value') and hasattr(self.config, 'final_gp_value') and hasattr(self.config, 'gp_step'):
            gp_str = f"{self.config.group_penalty_type}: {self.config.initial_gp_value} to {self.config.final_gp_value} (step: {self.config.gp_step})"
            self.add_summary_item(layout, "Group Penalty:", gp_str)
        
        # Logspace and other settings
        if hasattr(self.config, 'use_logspace'):
            logspace_str = f"{'Yes' if self.config.use_logspace else 'No'}"
            if self.config.use_logspace and hasattr(self.config, 'num_log_points'):
                logspace_str += f" ({self.config.num_log_points} points)"
            self.add_summary_item(layout, "Use Logspace:", logspace_str)
        
        # Phenotype names
        if hasattr(self.config, 'pheno_name1') and hasattr(self.config, 'pheno_name2'):
            self.add_summary_item(layout, "Phenotype Names:", f"{self.config.pheno_name1} vs {self.config.pheno_name2}")
        
        # Deletion canceler options
        if hasattr(self.config, 'nix_full_deletions') or hasattr(self.config, 'cancel_only_partner') or hasattr(self.config, 'min_pairs'):
            del_opts = []
            if hasattr(self.config, 'nix_full_deletions') and self.config.nix_full_deletions:
                del_opts.append("Exclude full deletions")
            if hasattr(self.config, 'cancel_only_partner') and not self.config.cancel_only_partner:
                del_opts.append("Cancel all partners")
            if hasattr(self.config, 'min_pairs') and self.config.min_pairs > 1:
                del_opts.append(f"Min pairs: {self.config.min_pairs}")
            
            if del_opts:
                self.add_summary_item(layout, "Deletion Options:", ", ".join(del_opts))
        
        # Output options
        if hasattr(self.config, 'output_dir') and self.config.output_dir:
            self.add_summary_item(layout, "Output Directory:", self.config.output_dir)
        
        if hasattr(self.config, 'output_file_base_name') and self.config.output_file_base_name:
            self.add_summary_item(layout, "Output Base Name:", self.config.output_file_base_name)
        
        # Output toggles
        if hasattr(self.config, 'keep_raw_output') or hasattr(self.config, 'show_selected_sites'):
            toggles = []
            if hasattr(self.config, 'keep_raw_output') and self.config.keep_raw_output:
                toggles.append("Keep raw output")
            if hasattr(self.config, 'show_selected_sites') and self.config.show_selected_sites:
                toggles.append("Show selected sites")
            
            if toggles:
                self.add_summary_item(layout, "Output Options:", ", ".join(toggles))
        
        # Plot options
        if hasattr(self.config, 'make_sps_plot') or hasattr(self.config, 'make_sps_kde_plot'):
            plots = []
            if hasattr(self.config, 'make_sps_plot') and self.config.make_sps_plot:
                plots.append("SPS Plot")
            if hasattr(self.config, 'make_sps_kde_plot') and self.config.make_sps_kde_plot:
                plots.append("SPS KDE Plot")
            
            if plots:
                self.add_summary_item(layout, "Generate Plots:", ", ".join(plots))
        
        # Multi-matrix options
        multi_opts = []
        if hasattr(self.config, 'top_rank_frac') and self.config.top_rank_frac != 0.01:
            multi_opts.append(f"Top rank frac: {self.config.top_rank_frac}")
        if hasattr(self.config, 'response_dir') and self.config.response_dir:
            multi_opts.append(f"Response dir: {self.config.response_dir}")
        
        # Null model options
        if hasattr(self.config, 'make_null_models') and self.config.make_null_models:
            multi_opts.append("Generate null models")
        if hasattr(self.config, 'make_pair_randomized_null_models') and self.config.make_pair_randomized_null_models:
            rand_str = "Generate pair-randomized nulls"
            if hasattr(self.config, 'num_randomized_alignments') and self.config.num_randomized_alignments != 10:
                rand_str += f" ({self.config.num_randomized_alignments} randomizations)"
            multi_opts.append(rand_str)
        
        if multi_opts:
            self.add_summary_item(layout, "Multi-matrix Options:", "; ".join(multi_opts))
    
    def add_summary_item(self, layout, label, value):
        """Add a single item to the configuration summary."""
        label_widget = QLabel(label)
        value_widget = QLineEdit(value)
        value_widget.setReadOnly(True)
        value_widget.setStyleSheet("background: #f8f8f8; border: 1px solid #ddd;")
        layout.addRow(label_widget, value_widget)
    
    def clear_layout(self, layout):
        """Clear all widgets from a layout."""
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    def run_analysis(self):
        """Run the ESL-PSC analysis."""
        # Create a new RunPage instance and show it
        run_page = RunPage(self.config)
        self.wizard().addPage(run_page)
        self.wizard().setCurrentPage(run_page)


class RunPage(BaseWizardPage):
    """Page for running the analysis and displaying progress."""
    
    def __init__(self, config):
        """Initialize the run page."""
        super().__init__("Run Analysis")
        self.config = config
        self.setSubTitle("Run the analysis and monitor progress.")
        
        # Thread pool for running the worker
        self.thread_pool = QThreadPool()
        self.worker = None
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
        # Command display section
        cmd_group = QGroupBox("ESL-PSC Command")
        cmd_layout = QVBoxLayout()
        
        # Command display with monospace font for better readability
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        self.cmd_display.setFont(QFont("Courier", 10))
        self.cmd_display.setPlaceholderText("The ESL-PSC command will be executed...")
        self.cmd_display.setMinimumHeight(100)
        
        # Copy button
        self.copy_btn = QPushButton("Copy to Clipboard")
        self.copy_btn.clicked.connect(self.copy_command_to_clipboard)
        
        # Add widgets to command layout
        cmd_layout.addWidget(QLabel("Command being executed:"))
        cmd_layout.addWidget(self.cmd_display)
        cmd_layout.addWidget(self.copy_btn, 0, Qt.AlignmentFlag.AlignRight)
        cmd_group.setLayout(cmd_layout)
        container_layout.addWidget(cmd_group)
        
        # Output display
        output_group = QGroupBox("Analysis Output")
        output_layout = QVBoxLayout()
        
        self.output_display = QTextEdit()
        self.output_display.setReadOnly(True)
        self.output_display.setFont(QFont("Courier", 10))
        self.output_display.setPlaceholderText("Analysis output will appear here...")
        output_layout.addWidget(self.output_display)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setTextVisible(True)
        output_layout.addWidget(QLabel("Progress:"))
        output_layout.addWidget(self.progress_bar)
        
        output_group.setLayout(output_layout)
        container_layout.addWidget(output_group, 1)  # Take remaining space
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.go_back)
        
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_analysis)
        
        btn_layout.addWidget(self.back_btn)
        btn_layout.addStretch()
        btn_layout.addWidget(self.run_btn)
        btn_layout.addWidget(self.stop_btn)
        
        container_layout.addLayout(btn_layout)
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)
    
    def copy_command_to_clipboard(self):
        """Copy the current command to the system clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.cmd_display.toPlainText())
    
    def go_back(self):
        """Go back to the previous page."""
        self.wizard().back()
    
    def on_enter(self):
        """Update the command display when the page is shown."""
        # Build the command (same as in OutputPage)
        cmd_parts = ["esl-psc"]
        
        # Required parameters
        if hasattr(self.config, 'alignment_dir') and self.config.alignment_dir:
            cmd_parts.append(f'--alignment-dir "{self.config.alignment_dir}"')
        if hasattr(self.config, 'species_groups_file') and self.config.species_groups_file:
            cmd_parts.append(f'--species-groups "{self.config.species_groups_file}"')
        
        # Optional parameters
        if hasattr(self.config, 'species_phenotypes_file') and self.config.species_phenotypes_file:
            cmd_parts.append(f'--species-phenotypes "{self.config.species_phenotypes_file}"')
        if hasattr(self.config, 'prediction_alignments_dir') and self.config.prediction_alignments_dir:
            cmd_parts.append(f'--prediction-alignments "{self.config.prediction_alignments_dir}"')
        if hasattr(self.config, 'limited_genes_file') and self.config.limited_genes_file:
            cmd_parts.append(f'--limited-genes "{self.config.limited_genes_file}"')
        if hasattr(self.config, 'output_dir') and self.config.output_dir:
            cmd_parts.append(f'--output-dir "{self.config.output_dir}"')
        
        # Join command parts with spaces
        cmd_str = ' \
  '.join(cmd_parts)
        
        # Display the command
        self.cmd_display.setPlainText(cmd_str)
        
        # Clear previous output
        self.output_display.clear()
        
        # Reset progress
        self.progress_bar.setValue(0)
        
        # Enable/disable buttons
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
    
    def run_analysis(self):
        """Run the ESL-PSC analysis."""
        try:
            cmd = self.config.get_command_string()
            
            # Update UI
            self.run_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            self.cmd_display.clear()
            self.cmd_display.append(f"$ {cmd}\n")
            
            # Create and configure worker
            self.worker = ESLWorker(cmd)
            self.worker.signals.output.connect(self.append_output)
            self.worker.signals.error.connect(self.append_error)
            self.worker.signals.finished.connect(self.analysis_finished)
            self.worker.signals.progress.connect(self.update_progress)
            
            # Start worker
            self.thread_pool.start(self.worker)
            
        except Exception as e:
            self.append_error(f"Error starting analysis: {str(e)}")
            self.analysis_finished(1)
    
    def stop_analysis(self):
        """Stop the running analysis."""
        if self.worker and self.worker.is_running:
            self.worker.stop()
            self.append_output("\nAnalysis stopped by user.")
            self.analysis_finished(-1)
    
    def append_output(self, text):
        """Append text to the output display."""
        self.cmd_display.append(text)
        self.cmd_display.verticalScrollBar().setValue(
            self.cmd_display.verticalScrollBar().maximum()
        )
    
    def append_error(self, text):
        """Append error text to the output display."""
        self.cmd_display.append(f"<span style='color:red'>{text}</span>")
        self.cmd_display.verticalScrollBar().setValue(
            self.cmd_display.verticalScrollBar().maximum()
        )
    
    def update_progress(self, value):
        """Update the progress bar."""
        self.progress_bar.setValue(value)
    
    def analysis_finished(self, exit_code):
        """Handle analysis completion."""
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        
        if exit_code == 0:
            self.append_output("\nAnalysis completed successfully!")
        elif exit_code == -1:
            self.append_output("\nAnalysis was stopped.")
        else:
            self.append_output(f"\nAnalysis failed with exit code {exit_code}.")
