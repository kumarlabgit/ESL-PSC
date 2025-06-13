#!/usr/bin/env python3
"""
ESL-PSC GUI Wizard

A graphical interface for running ESL-PSC analyses with an intuitive wizard.
"""
import sys
import os
from PyQt6.QtWidgets import QApplication, QMessageBox
from PyQt6.QtCore import Qt
from gui.ui.main_window import MainWindow

def main():
    """Main entry point for the ESL-PSC GUI application."""
    print("Starting ESL-PSC GUI...")
    print(f"Python executable: {sys.executable}")
    print(f"Working directory: {os.getcwd()}")
    
    try:
        # Try to enable high DPI scaling if available
        try:
            QApplication.setAttribute(Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
            QApplication.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, True)
            print("Enabled high DPI scaling")
        except (AttributeError, NameError) as e:
            print(f"Could not enable high DPI scaling (not supported in this Qt version): {e}")
        
        # Create application instance
        print("Creating QApplication...")
        app = QApplication(sys.argv)
        print("QApplication created")
        print("QApplication created")
        
        # Set application metadata
        app.setApplicationName("ESL-PSC Wizard")
        app.setApplicationVersion("0.1.0")
        app.setOrganizationName("ESL-PSC")
        
        # Create and show main window
        print("Creating main window...")
        window = MainWindow()
        print("Main window created")
        
        # Show the window
        window.show()
        print("Main window shown")
        
        # Show a test message box to verify GUI is working
        QMessageBox.information(
            window, 
            "ESL-PSC GUI", 
            "Application started successfully!\n\n"
            "If you can see this message, the GUI is working."
        )
        
        print("Starting event loop...")
        # Start the event loop
        sys.exit(app.exec())
        
    except Exception as e:
        print(f"FATAL ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Try to show error in a message box if possible
        try:
            error_msg = f"Fatal error: {str(e)}\n\n{traceback.format_exc()}"
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Critical)
            msg.setText("Fatal Error")
            msg.setInformativeText(str(e))
            msg.setDetailedText(traceback.format_exc())
            msg.setWindowTitle("ESL-PSC GUI Error")
            msg.exec()
        except:
            pass  # If we can't show the error in a GUI, we've already printed it
            
        sys.exit(1)

if __name__ == "__main__":
    main()
