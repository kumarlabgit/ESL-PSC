"""
Main window for the ESL-PSC Wizard application.
"""
import os
from PyQt6.QtWidgets import (
    QMainWindow, QWizard, QWizardPage, QVBoxLayout, QLabel, QWidget,
    QGroupBox, QFormLayout, QScrollArea, QSizePolicy, QSpacerItem,
    QDoubleSpinBox, QSpinBox, QComboBox, QCheckBox, QMessageBox,
    QTextEdit, QLineEdit, QPushButton, QProgressBar, QHBoxLayout, QApplication
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
            self.setOption(QWizard.WizardOption.HaveFinishButtonOnEarlyPages, True)
            self.setOption(QWizard.WizardOption.NoBackButtonOnLastPage, True)
            self.setOption(QWizard.WizardOption.NoCancelButton, False)
            self.setMinimumSize(900, 650)
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
        req_layout = QFormLayout()
        
        # Alignment directory
        self.alignment_dir = FileSelector(
            "Alignment Directory:", 'directory',
            default_path=os.getcwd()
        )
        self.alignment_dir.path_changed.connect(
            lambda p: setattr(self.config, 'alignment_dir', p)
        )
        req_layout.addRow(self.alignment_dir)
        
        # Species groups file
        self.species_groups = FileSelector(
            "Species Groups File:", 'file',
            default_path=os.getcwd()
        )
        self.species_groups.path_changed.connect(
            lambda p: setattr(self.config, 'species_groups_file', p)
        )
        req_layout.addRow(self.species_groups)
        
        req_group.setLayout(req_layout)
        container_layout.addWidget(req_group)
        
        # Optional inputs group
        opt_group = QGroupBox("Optional Inputs")
        opt_layout = QFormLayout()
        
        # Species phenotypes file
        self.species_phenotypes = FileSelector(
            "Species Phenotypes File:", 'file',
            default_path=os.getcwd()
        )
        self.species_phenotypes.path_changed.connect(
            lambda p: setattr(self.config, 'species_phenotypes_file', p)
        )
        opt_layout.addRow(self.species_phenotypes)
        
        # Prediction alignments directory
        self.prediction_alignments = FileSelector(
            "Prediction Alignments Directory:", 'directory',
            default_path=os.getcwd()
        )
        self.prediction_alignments.path_changed.connect(
            lambda p: setattr(self.config, 'prediction_alignments_dir', p)
        )
        opt_layout.addRow(self.prediction_alignments)
        
        # Limited genes file
        self.limited_genes = FileSelector(
            "Limited Genes File:", 'file',
            default_path=os.getcwd()
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
    
    def __init__(self, config):
        """Initialize the parameters page."""
        super().__init__("Analysis Parameters")
        self.config = config
        self.setSubTitle("Configure the analysis parameters.")
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
        # Regularization parameters group
        reg_group = QGroupBox("Regularization Parameters")
        reg_layout = QFormLayout()
        
        # Sparsity parameters
        self.lambda1 = QDoubleSpinBox()
        self.lambda1.setRange(0.0, 1.0)
        self.lambda1.setSingleStep(0.01)
        self.lambda1.setValue(0.1)
        self.lambda1.valueChanged.connect(
            lambda v: setattr(self.config, 'lambda1', v)
        )
        reg_layout.addRow("Lambda 1:", self.lambda1)
        
        self.lambda2 = QDoubleSpinBox()
        self.lambda2.setRange(0.0, 1.0)
        self.lambda2.setSingleStep(0.01)
        self.lambda2.setValue(0.1)
        self.lambda2.valueChanged.connect(
            lambda v: setattr(self.config, 'lambda2', v)
        )
        reg_layout.addRow("Lambda 2:", self.lambda2)
        
        # Group penalties
        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(0.0, 1.0)
        self.alpha.setSingleStep(0.1)
        self.alpha.setValue(0.5)
        self.alpha.valueChanged.connect(
            lambda v: setattr(self.config, 'alpha', v)
        )
        reg_layout.addRow("Alpha:", self.alpha)
        
        reg_group.setLayout(reg_layout)
        container_layout.addWidget(reg_group)
        
        # Filtering options group
        filter_group = QGroupBox("Filtering Options")
        filter_layout = QFormLayout()
        
        # Deletion handling
        self.deletion_handling = QComboBox()
        self.deletion_handling.addItems(["none", "gap-canceled", "complete"])
        self.deletion_handling.setCurrentText("gap-canceled")
        self.deletion_handling.currentTextChanged.connect(
            lambda t: setattr(self.config, 'deletion_handling', t)
        )
        filter_layout.addRow("Deletion Handling:", self.deletion_handling)
        
        # Min species per gene
        self.min_species = QSpinBox()
        self.min_species.setRange(1, 1000)
        self.min_species.setValue(4)
        self.min_species.valueChanged.connect(
            lambda v: setattr(self.config, 'min_species_per_gene', v)
        )
        filter_layout.addRow("Min Species per Gene:", self.min_species)
        
        filter_group.setLayout(filter_layout)
        container_layout.addWidget(filter_group)
        
        # Output options group
        output_group = QGroupBox("Output Options")
        output_layout = QFormLayout()
        
        # Output directory
        self.output_dir = FileSelector(
            "Output Directory:", 'directory',
            default_path=os.path.join(os.getcwd(), 'esl_psc_output')
        )
        self.output_dir.path_changed.connect(
            lambda p: setattr(self.config, 'output_dir', p)
        )
        output_layout.addRow(self.output_dir)
        
        # Plot format
        self.plot_format = QComboBox()
        self.plot_format.addItems(["pdf", "png", "svg"])
        self.plot_format.setCurrentText("pdf")
        self.plot_format.currentTextChanged.connect(
            lambda t: setattr(self.config, 'plot_format', t)
        )
        output_layout.addRow("Plot Format:", self.plot_format)
        
        output_group.setLayout(output_layout)
        container_layout.addWidget(output_group)
        
        # Add stretch to push everything to the top
        container_layout.addStretch()
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)


class OutputPage(BaseWizardPage):
    """Page for reviewing configuration before running the analysis."""
    
    def __init__(self, config):
        """Initialize the output page."""
        super().__init__("Review Configuration")
        self.config = config
        self.setSubTitle("Review your configuration before running the analysis.")
        
        # Command preview
        self.cmd_preview = QTextEdit()
        self.cmd_preview.setReadOnly(True)
        self.cmd_preview.setFont(QFont("Courier"))
        self.cmd_preview.setPlaceholderText("Configuration summary will appear here...")
        
        # Add to layout - using the layout from BaseWizardPage
        self.layout().addWidget(QLabel("ESL-PSC Command:"))
        self.layout().addWidget(self.cmd_preview)
        
        # Add stretch to push everything to the top
        self.layout().addStretch()
    
    def on_enter(self):
        """Update the command preview when the page is shown."""
        try:
            cmd = self.config.get_command_string()
            self.cmd_preview.setPlainText(cmd)
        except Exception as e:
            self.cmd_preview.setPlainText(f"Error generating command: {str(e)}")


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
        
        # Command display
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        self.cmd_display.setFont(QFont("Courier"))
        self.cmd_display.setPlaceholderText("Command output will appear here...")
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setTextVisible(True)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_analysis)
        
        btn_layout.addWidget(self.run_btn)
        btn_layout.addWidget(self.stop_btn)
        
        # Add widgets to layout - using the layout from BaseWizardPage
        self.layout().addWidget(QLabel("Analysis Output:"))
        self.layout().addWidget(self.cmd_display)
        self.layout().addWidget(QLabel("Progress:"))
        self.layout().addWidget(self.progress_bar)
        self.layout().addLayout(btn_layout)
    
    def on_enter(self):
        """Update the command display when the page is shown."""
        try:
            cmd = self.config.get_command_string()
            self.cmd_display.setPlainText(f"$ {cmd}\n\n")
        except Exception as e:
            self.cmd_display.setPlainText(f"Error: {str(e)}")
    
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
