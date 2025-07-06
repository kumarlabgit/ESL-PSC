"""
File and directory selection widgets for the ESL-PSC GUI.
"""
import os
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog
)
from PySide6.QtCore import Signal

class FileSelector(QWidget):
    """A widget for selecting a file or directory with a browse button and optional description."""
    
    path_changed = Signal(str)  # Signal emitted when path changes
    # Class-level attribute tracking the most recently selected directory across all widgets
    last_selected_dir = os.getcwd()
    
    def __init__(self, label, mode='file', default_path='', description='', parent=None):
        """
        Initialize the file selector.
        
        Args:
            label: The label to display next to the path field
            mode: Either 'file' or 'directory'
            default_path: The default path to display
            description: Optional help text to display below the input
            parent: The parent widget
        """
        super().__init__(parent)
        self.mode = mode
        self.default_path = default_path or FileSelector.last_selected_dir
        
        # Create main layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(2)  # Reduced spacing between widgets
        
        # Create top row layout (label + input + button)
        top_row = QHBoxLayout()
        top_row.setContentsMargins(0, 0, 0, 0)
        
        # Add label
        top_row.addWidget(QLabel(label))
        
        # Add path display
        self.path_edit = QLineEdit()
        self.path_edit.setReadOnly(False)
        self.path_edit.setPlaceholderText(f"Select {mode}...")
        # Allow manual editing and propagate changes when editing finishes
        self.path_edit.editingFinished.connect(self._on_edit_finished)
        top_row.addWidget(self.path_edit, 1)  # Stretch factor 1
        
        # Add browse button
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse)
        top_row.addWidget(self.browse_btn)
        
        main_layout.addLayout(top_row)
        
        # Add description label if provided
        if description:
            self.desc_label = QLabel(description)
            self.desc_label.setWordWrap(True)
            self.desc_label.setStyleSheet("color: #666666; font-size: 10pt;")
            self.desc_label.setContentsMargins(5, 0, 0, 5)  # Left indent for description
            main_layout.addWidget(self.desc_label)
    
    def browse(self):
        """Open a file/directory dialog to select the path.

        The dialog starts in the directory that was most recently selected in any
        FileSelector instance.  This greatly speeds up workflows that require
        choosing multiple related files.
        """
        start_dir = FileSelector.last_selected_dir or self.default_path

        if self.mode == 'file':
            path, _ = QFileDialog.getOpenFileName(
                self, f"Select {self.mode}", start_dir
            )
        else:  # directory
            path = QFileDialog.getExistingDirectory(
                self, f"Select {self.mode}", start_dir
            )

        if path:
            self.set_path(path)
    
    def get_path(self):
        """Get the currently selected path."""
        return self.path_edit.text()
    
    def set_path(self, path):
        """Set the current path and emit the changed signal.

        An empty string should clear the field instead of displaying '.'.
        """
        if not path:                             # handle '' or None cleanly
            self.path_edit.clear()
            self.path_changed.emit("")
            return

        # Normal non-empty path
        path = os.path.normpath(path)
        self.path_edit.setText(path)
        self.path_changed.emit(path)

        # Update default path for both this instance and all other FileSelector instances
        if os.path.exists(path):
            # Always use the PARENT directory of the selected item as the next
            # default location.  This applies to both files and directories so
            # that selecting a folder does not cause the dialog to reopen
            # *inside* that folder.
            next_dir = os.path.dirname(path)
            # Edge-case: os.path.dirname('/') returns '/', keep it.
            self.default_path = next_dir or path
            FileSelector.last_selected_dir = self.default_path
    
    def _on_edit_finished(self):
        """Synchronise manual line-edit changes with the widget state."""
        self.set_path(self.path_edit.text())

    def clear(self):
        """Clear the current selection."""
        self.path_edit.clear()
        self.path_changed.emit("")
