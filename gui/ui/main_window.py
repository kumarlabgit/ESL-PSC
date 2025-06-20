"""
Main window for the ESL-PSC Wizard application.
"""
import os
from PyQt6.QtWidgets import (
    QMainWindow, QWizard, QVBoxLayout, QWidget, QPushButton, QApplication,
    QFileDialog, QMessageBox
)

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
        # Width of the Next button used as an invisible placeholder on the last page
        self._next_btn_placeholder_width = None
        
        try:
            # Set window properties first
            self.setObjectName("ESLWizard")
            self.setWindowTitle("ESL-PSC Analysis Wizard")
            self.setWizardStyle(QWizard.WizardStyle.ModernStyle)
            self.setOption(QWizard.WizardOption.HaveHelpButton, False)
            # Hide the Back button completely on the first page for consistency
            self.setOption(QWizard.WizardOption.NoBackButtonOnStartPage, True)
            # Keep the (disabled) Next button present on the last page so that the
            # Back button does not shift position when the user navigates there.
            # We will later disable and blank the button, but leaving it in the
            # layout preserves the spacing that users rely on.
            self.setOption(QWizard.WizardOption.HaveNextButtonOnLastPage, True)
            self.setOption(QWizard.WizardOption.HaveFinishButtonOnEarlyPages, False)  # No grayed-out finish button
            self.setOption(QWizard.WizardOption.NoBackButtonOnLastPage, False)
            self.setOption(QWizard.WizardOption.NoCancelButton, False)
            self.setMinimumSize(800, 700)  # Reduced minimum width to align with MainWindow
            print("ESLWizard: Window properties set")

            # Create custom Save/Load buttons before setting cancel button text
            save_btn = QPushButton("Save Config")
            load_btn = QPushButton("Load Config")
            save_btn.clicked.connect(self.save_config)
            load_btn.clicked.connect(self.load_config)
            # Assign to custom button slots and define layout (left-aligned)
            self.setButton(QWizard.WizardButton.CustomButton1, save_btn)
            self.setButton(QWizard.WizardButton.CustomButton2, load_btn)
            self.setButtonLayout([
                QWizard.WizardButton.CustomButton1,
                QWizard.WizardButton.CustomButton2,
                QWizard.WizardButton.Stretch,
                QWizard.WizardButton.BackButton,
                QWizard.WizardButton.NextButton,
                QWizard.WizardButton.CancelButton,
            ])

            self.setButtonText(QWizard.WizardButton.CancelButton, "Quit")
            quit_btn = self.button(QWizard.WizardButton.CancelButton)
            try:
                quit_btn.clicked.disconnect()
            except TypeError:
                pass
            quit_btn.clicked.connect(self.confirm_quit)
            
            # Apply stylesheet once now (it will be refreshed on any palette change)
            self.apply_stylesheet()
            # Enable automatic dark / light switching
            self._connect_color_scheme_updates()
            
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
    
    # ──────────────────────────────────────────────────────────────────────
    #   Dynamic dark / light‑mode support
    # ──────────────────────────────────────────────────────────────────────
    def apply_stylesheet(self):
        """(Re)apply the wizard‑wide stylesheet using the current QPalette."""
        # Determine border color for inputs based on theme, as palette(mid) is unreliable
        # --- Robust dark‑mode detection (Qt 6.5+ palette‑based) -------------
        pal = QApplication.palette()
        window_lum = pal.color(QPalette.ColorRole.Window).lightness()
        text_lum   = pal.color(QPalette.ColorRole.WindowText).lightness()
        is_dark_mode = window_lum < text_lum          # darker background → dark theme
        border_color = "rgba(255, 255, 255, 180)" if is_dark_mode else "rgba(0, 0, 0, 80)"

        # DEBUG: print once every time the stylesheet is (re)applied
        print(f"[apply_stylesheet] dark_mode={is_dark_mode}  "
              f"window_lum={window_lum}  text_lum={text_lum}  border={border_color}")

        # Build the stylesheet with a placeholder token then replace it with the
        # actual border color.  This avoids the heavy brace-escaping needed in
        # f-strings when embedding large chunks of QSS.
        css = """
            /* Use the current system palette for backgrounds */
            QWizard            {{ background: palette(window); }}
            QWizardPage        {{ background: palette(window); padding: 20px; }}

            /* Page labels/text */
            QWizardPage > QLabel {{
                font-size: 14px;
                margin-bottom: 10px;
                color: palette(text);
            }}

            /* Section group‑boxes */
            QWizardPage > QGroupBox {{
                font-weight: bold;
                border: 1px solid palette(mid);
                border-radius: 5px;
                margin-top: 2ex;
                padding: 10px;
                color: palette(text);
            }}
            QWizardPage > QGroupBox::title {{
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 3px 0 3px;
                color: palette(text);
            }}

            /* General input widgets – give them a visible border */
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit {
                background-color: palette(base);
                color: palette(text);
                border: 1px solid __BORDER__;
                border-radius: 3px;
                padding: 4px;
            }
            QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QTextEdit:focus, QComboBox:focus {{
                border: 1px solid palette(highlight);
            }}

            /* Read-only summary fields */
            QLineEdit[readOnly="true"] {{
                background-color: palette(window);
                color: palette(text);
                border: 1px solid __BORDER__;
            }}

            /* Use built-in plus/minus buttons for spinboxes */
            QSpinBox, QDoubleSpinBox {{
                qproperty-buttonSymbols: UpDownArrows;
            }}

            /* Command-summary group labels & fields */
            QGroupBox#configSummaryGroup > QLabel {{ color: palette(text); }}
            QGroupBox#configSummaryGroup QLineEdit {{
                background-color: palette(window);
                color: palette(text);
                border: 1px solid __BORDER__;
                padding: 4px;
            }}
        """
        # Substitute the placeholder token with the real colour
        css = css.replace("__BORDER__", border_color)

        # Apply to the entire application so every widget—present and future—
        # inherits the same QSS, and dark-mode borders cannot be overridden
        # later by Qt’s internal style.
        qApp = QApplication.instance()
        if qApp:
            qApp.setStyleSheet(css)
            # Force a full refresh of the style system
            qApp.style().polish(qApp)

    # --- Automatic palette-change support ----------------------------------
    # --- Automatic palette-change support ----------------------------------
    def _connect_color_scheme_updates(self):
        """
        Re-apply the palette-aware stylesheet whenever the OS theme flips.
        We listen to two independent signals so the code works on every
        platform / Qt build:
          1) QApplication.instance().paletteChanged  – emits on macOS, Wayland
          2) QStyleHints.colorSchemeChanged          – emits on X11, Windows
        """
        app = QApplication.instance()
        if app and hasattr(app, "paletteChanged"):
            app.paletteChanged.connect(self.apply_stylesheet)

        sh = QApplication.styleHints()
        if hasattr(sh, "colorSchemeChanged"):
            sh.colorSchemeChanged.connect(self.apply_stylesheet)

    # Fallback event-filter – only used on Qt < 6.3
    def eventFilter(self, obj, event):
        if event.type() == QEvent.Type.ApplicationPaletteChange:
            if not getattr(self, "_palette_guard", False):
                self._palette_guard = True          # prevent recursion
                self.apply_stylesheet()
                self._palette_guard = False
        return super().eventFilter(obj, event)

    # ──────────────────────────────────────────────────────────────────────────
    # Save / Load configuration helpers
    # ──────────────────────────────────────────────────────────────────────────
    def save_config(self):
        """Prompt user to save current configuration to a JSON file."""
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Configuration",
            os.getcwd(),
            "JSON Files (*.json)"
        )
        if path:
            # Ensure .json extension
            if not path.lower().endswith('.json'):
                path += '.json'
            try:
                with open(path, "w", encoding="utf-8") as fh:
                    json.dump(asdict(self.config), fh, indent=2)
            except Exception as exc:
                QMessageBox.critical(self, "Save Error", f"Could not save configuration:\n{exc}")

    def load_config(self):
        """Prompt user to load configuration from a JSON file and update UI."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Configuration",
            os.getcwd(),
            "JSON Files (*.json)"
        )
        if path:
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    data = json.load(fh)
                # Update attributes that exist on ESLConfig
                for key, val in data.items():
                    if hasattr(self.config, key):
                        setattr(self.config, key, val)
                self.refresh_pages_from_config()
            except Exception as exc:
                QMessageBox.critical(self, "Load Error", f"Could not load configuration:\n{exc}")

    def refresh_pages_from_config(self):
        """Notify pages to sync their widgets with the current config."""
        for page in [self.input_page, self.params_page, self.command_page]:
            if hasattr(page, "update_ui_from_config"):
                page.update_ui_from_config()
        # If currently on the command page, refresh its display immediately
        if self.currentPage() == self.command_page:
            self.command_page.on_enter()
    
    def on_current_id_changed(self, page_id):
        """Handle page changes in the wizard."""
        current_page = self.currentPage()
        if hasattr(current_page, 'on_enter'):
            current_page.on_enter()

        # Keep Back button from shifting when on last page.
        next_btn = self.button(QWizard.WizardButton.NextButton)
        if next_btn is None:
            return  # Safety – should not happen

        # Cache the natural width of a normal "Next" button the first time
        if self._next_btn_placeholder_width is None:
            self._next_btn_placeholder_width = next_btn.sizeHint().width()

        # If on the last page, disable and visually blank the Next button so it
        # takes up space but is not interactable or confusing to users.
        if self.currentId() == self.pageIds()[-1]:
            next_btn.setEnabled(False)
            next_btn.setText("")
            next_btn.setFixedWidth(self._next_btn_placeholder_width)
            # Remove borders to make it appear invisible while still reserving space
            next_btn.setStyleSheet("border: none; background: transparent;")
        else:
            # Restore the normal appearance/behaviour on non-last pages
            next_btn.setEnabled(True)
            next_btn.setFixedWidth(self._next_btn_placeholder_width)
            next_btn.setStyleSheet("")
            if next_btn.text() == "":
                next_btn.setText("Next")

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

            # Depending on which input type is selected, require the correct path
            if self.input_page.use_species_groups.isChecked():
                if not self.config.species_groups_file:
                    QMessageBox.warning(self, "Missing Required Field",
                                        "Please select a species groups file.")
                    return False
            else:
                if not self.config.response_dir:
                    QMessageBox.warning(self, "Missing Required Field",
                                        "Please select a response matrix directory.")
                    return False

        return super().validateCurrentPage()


# --- pulled-out page classes live in gui.ui.pages ---
from gui.ui.pages.input_page import InputPage
from gui.ui.pages.parameters_page import ParametersPage
from gui.ui.pages.command_page import CommandPage
from gui.ui.pages.run_page import RunPage
