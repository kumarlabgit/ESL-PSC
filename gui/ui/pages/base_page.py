"""Common functionality shared by all wizard pages."""
from PySide6.QtWidgets import QWizardPage, QVBoxLayout

class BaseWizardPage(QWizardPage):
    """Base class for wizard pages with common functionality."""

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setTitle(title)
        # Remove the default margins that Qt adds around QWizardPage so the
        # page’s contents can span the full available width.
        self.setContentsMargins(0, 0, 0, 0)
        # Only set a layout once
        if not self.layout():
            self.setLayout(QVBoxLayout())

        # Eliminate default page margins so page contents can span full width.
        # Keep zero margins on all sides – individual child layouts/group-boxes
        # already handle their own internal spacing.
        self.layout().setContentsMargins(0, 0, 0, 0)

    # Sub-classes can override this hook
    def on_enter(self):
        """Called automatically when the page becomes current."""
        pass
