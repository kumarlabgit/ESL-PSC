"""
File and directory selection widgets for the ESL-PSC GUI.
"""
import os
from pathlib import Path
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QMessageBox
)
from PyQt6.QtCore import Qt, pyqtSignal

class FileSelector(QWidget):
    """A widget for selecting a file or directory with a browse button."""
    
    path_changed = pyqtSignal(str)  # Signal emitted when path changes
    
    def __init__(self, label, mode='file', default_path='', parent=None):
        """
        Initialize the file selector.
        
        Args:
            label: The label to display next to the path field
            mode: Either 'file' or 'directory'
            default_path: The default path to display
            parent: The parent widget
        """
        super().__init__(parent)
        self.mode = mode
        self.default_path = default_path or os.getcwd()
        
        # Create layout
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Add label
        layout.addWidget(QLabel(label))
        
        # Add path display
        self.path_edit = QLineEdit()
        self.path_edit.setReadOnly(True)
        self.path_edit.setPlaceholderText(f"Select {mode}...")
        layout.addWidget(self.path_edit, 1)  # Stretch factor 1
        
        # Add browse button
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse)
        layout.addWidget(self.browse_btn)
    
    def browse(self):
        """Open a file/directory dialog to select the path."""
        if self.mode == 'file':
            path, _ = QFileDialog.getOpenFileName(
                self, f"Select {self.mode}", self.default_path
            )
        else:  # directory
            path = QFileDialog.getExistingDirectory(
                self, f"Select {self.mode}", self.default_path
            )
        
        if path:
            self.set_path(path)
    
    def get_path(self):
        """Get the currently selected path."""
        return self.path_edit.text()
    
    def set_path(self, path):
        """Set the current path and emit the changed signal."""
        path = os.path.normpath(path)
        self.path_edit.setText(path)
        self.path_changed.emit(path)
        
        # Update default path to the parent directory
        if os.path.exists(path):
            if os.path.isfile(path):
                self.default_path = os.path.dirname(path)
            else:
                self.default_path = path
    
    def clear(self):
        """Clear the current selection."""
        self.path_edit.clear()
        self.path_changed.emit("")
