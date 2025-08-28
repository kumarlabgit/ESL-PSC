from __future__ import annotations

from PySide6.QtWidgets import QMessageBox


class AutoSelectOptionsDialog(QMessageBox):
    """Simple dialog to choose auto-select behavior.

    Presents options in this order: Longest Sequence → Random → Default → Cancel.
    Stores the chosen option as `self.choice` in {"default","longest","random",None}.
    """

    def __init__(self, allow_longest: bool, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Auto Select Options")
        self.setText("Choose how to resolve equally valid siblings")

        # Buttons in a fixed order using ActionRole to preserve ordering cross-platform
        self.longest_btn = self.addButton(
            "Longest Sequence", QMessageBox.ButtonRole.ActionRole
        )
        if not allow_longest:
            self.longest_btn.setEnabled(False)
            self.longest_btn.setToolTip(
                "Set an alignments directory in the input page, to enable this option"
            )

        self.random_btn = self.addButton(
            "Random", QMessageBox.ButtonRole.ActionRole
        )

        self.default_btn = self.addButton(
            "Default", QMessageBox.ButtonRole.ActionRole
        )
        # Make the Default button appear as the default (blue) accept button.
        self.setDefaultButton(self.default_btn)

        # Cancel button
        self.cancel_btn = self.addButton(
            "Cancel", QMessageBox.ButtonRole.RejectRole
        )

        self.choice: str | None = None

    def exec(self) -> int:
        result = super().exec()
        clicked = self.clickedButton()
        if clicked == self.default_btn:
            self.choice = "default"
        elif clicked == self.longest_btn:
            self.choice = "longest"
        elif clicked == self.random_btn:
            self.choice = "random"
        else:
            self.choice = None
        return result

    def reject(self) -> None:
        """Handle programmatic rejection (e.g., Esc)."""
        self.choice = None
        super().reject()

    def closeEvent(self, event):
        """Treat window close as cancel."""
        self.choice = None
        super().closeEvent(event)
