"""
Main window for the ESL-PSC Wizard application.
"""
import os
from PyQt6.QtWidgets import (
    QMainWindow, QWizard, QWizardPage, QVBoxLayout, QLabel, QWidget,
    QGroupBox, QFormLayout, QScrollArea, QSizePolicy, QSpacerItem, QFrame, QStackedWidget,
    QDoubleSpinBox, QSpinBox, QComboBox, QCheckBox, QMessageBox, QButtonGroup,
    QTextEdit, QLineEdit, QPushButton, QProgressBar, QHBoxLayout, QApplication, QRadioButton,
    QFileDialog, QAbstractSpinBox
)
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt, pyqtSignal, QThreadPool
from PyQt6.QtGui import QFont, QTextCursor

from gui.core.config import ESLConfig

class MainWindow(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        """Initialize the main window and set up the wizard."""
        super().__init__()
        print("MainWindow: Initializing...")
        
        try:
            # Set window properties
            self.setWindowTitle("ESL-PSC Wizard")
            self.setMinimumSize(800, 900)  # Narrower and taller window
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
            self.setOption(QWizard.WizardOption.NoBackButtonOnLastPage, False)
            self.setOption(QWizard.WizardOption.NoCancelButton, False)
            self.setMinimumSize(1000, 700)  # Slightly larger minimum size
            print("ESLWizard: Window properties set")

            self.setButtonText(QWizard.WizardButton.CancelButton, "Quit")
            quit_btn = self.button(QWizard.WizardButton.CancelButton)
            try:
                quit_btn.clicked.disconnect()
            except TypeError:
                pass
            quit_btn.clicked.connect(self.confirm_quit)
            
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
            self.command_page = CommandPage(self.config)
            self.run_page = RunPage(self.config)
            
            # Add pages
            print("ESLWizard: Adding pages...")
            self.addPage(self.input_page)
            print("ESLWizard: Added InputPage")
            self.addPage(self.params_page)
            print("ESLWizard: Added ParametersPage")
            self.addPage(self.command_page)
            print("ESLWizard: Added CommandPage")
            self.addPage(self.run_page)
            print("ESLWizard: Added RunPage")

            # Hide the default “Finish” button – we use only Back / Next / Quit
            self.button(QWizard.WizardButton.FinishButton).hide()
            
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

    def confirm_quit(self):
        """
        Ask the user to confirm quitting. Called when the Quit button is pressed.
        """
        reply = QMessageBox.question(
            self,
            "Quit ESL-PSC Wizard",
            "Are you sure you want to quit?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            QApplication.instance().quit()

    # When the user presses Esc or some code calls reject(), show the same dialog
    def reject(self):
        self.confirm_quit()
    
    def validateCurrentPage(self):
        """Validate the current page before allowing the user to proceed."""
        current_page = self.currentPage()
        
        # Check required fields on input page
        if current_page == self.input_page:
            if not self.config.alignments_dir:
                QMessageBox.warning(self, "Missing Required Field", 
                                  "Please select an alignment directory.")
                return False
            if not self.config.species_groups_file:
                QMessageBox.warning(self, "Missing Required Field",
                                  "Please select a species groups file.")
                return False
        
        return super().validateCurrentPage()


# --- pulled-out page classes live in gui.ui.pages ---
from gui.ui.pages.base_page import BaseWizardPage
from gui.ui.pages.input_page import InputPage
from gui.ui.pages.parameters_page import ParametersPage
from gui.ui.pages.command_page import CommandPage
from gui.ui.pages.run_page import RunPage
