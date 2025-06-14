"""Common functionality shared by all wizard pages."""
from PyQt6.QtWidgets import QWizardPage, QVBoxLayout

class BaseWizardPage(QWizardPage):
    """Base class for wizard pages with common functionality."""

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setTitle(title)
        # Only set a layout once
        if not self.layout():
            self.setLayout(QVBoxLayout())

    # Sub-classes can override this hook
    def on_enter(self):
        """Called automatically when the page becomes current."""
        pass
