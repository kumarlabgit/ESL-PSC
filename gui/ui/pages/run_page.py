"""Run-analysis page of the ESL-PSC wizard (live output & progress)."""
from __future__ import annotations

from PyQt6.QtCore import QThreadPool, Qt
import os
from PyQt6.QtGui import QFontDatabase
from PyQt6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QPlainTextEdit, QPushButton,
    QLabel, QProgressBar, QHBoxLayout, QWizard, QMessageBox
)

from gui.core.worker import ESLWorker
from .base_page import BaseWizardPage
from gui.ui.widgets.results_display import (
    SpsPlotDialog, GeneRanksDialog, SelectedSitesDialog
)

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
        
        # Create scroll area and main container
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        container = QWidget()
        scroll.setWidget(container)
        container_layout = QVBoxLayout(container)
        
        # Terminal Output GroupBox
        cmd_group = QGroupBox("Analysis Terminal Output")
        cmd_layout = QVBoxLayout(cmd_group)
        
        self.cmd_display = QPlainTextEdit()
        self.cmd_display.setReadOnly(True)
        # Allow copying / selection on macOS and other platforms
        self.cmd_display.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.cmd_display.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse | Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        # Use a robust method to find the system's preferred monospace font
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setPointSize(10)
        self.cmd_display.setFont(font)
        # Placeholder text is supported by QPlainTextEdit starting Qt 6.2+
        if hasattr(self.cmd_display, "setPlaceholderText"):
            self.cmd_display.setPlaceholderText("Click 'Run Analysis' to start the ESL-PSC process...")
        cmd_layout.addWidget(self.cmd_display)
        container_layout.addWidget(cmd_group)

        # Result display buttons (hidden until analysis completes)
        self.results_layout = QHBoxLayout()
        self.sps_btn = QPushButton("Show SPS Plot")
        self.sps_btn.hide()
        self.sps_btn.clicked.connect(self.show_sps_plot)
        self.gene_btn = QPushButton("Show Top Gene Ranks")
        self.gene_btn.hide()
        self.gene_btn.clicked.connect(self.show_gene_ranks)
        self.sites_btn = QPushButton("Show Selected Sites")
        self.sites_btn.hide()
        self.sites_btn.clicked.connect(self.show_selected_sites)
        self.results_layout.addStretch()
        self.results_layout.addWidget(self.sps_btn)
        self.results_layout.addWidget(self.gene_btn)
        self.results_layout.addWidget(self.sites_btn)
        container_layout.addLayout(self.results_layout)

        # Progress Bars GroupBox
        progress_group = QGroupBox("Progress")
        progress_layout = QVBoxLayout(progress_group)

        # Overall Progress Bar
        progress_layout.addWidget(QLabel("Overall Progress (Species Combinations):"))
        self.overall_progress_bar = QProgressBar()
        self.overall_progress_bar.setRange(0, 100)
        self.overall_progress_bar.setTextVisible(True)
        progress_layout.addWidget(self.overall_progress_bar)

        # Step Status Label
        self.step_status_label = QLabel("Ready to run.")
        self.step_status_label.setStyleSheet("margin-top: 8px; font-style: italic;")
        progress_layout.addWidget(self.step_status_label)
        
        # Step Progress Bar
        self.step_progress_bar = QProgressBar()
        self.step_progress_bar.setRange(0, 100)
        self.step_progress_bar.setTextVisible(True)
        self.step_progress_bar.setFormat("Current Step: %p%")
        progress_layout.addWidget(self.step_progress_bar)

        container_layout.addWidget(progress_group)

        # Buttons
        btn_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_analysis)
        btn_layout.addStretch()
        btn_layout.addWidget(self.run_btn)
        btn_layout.addWidget(self.stop_btn)
        container_layout.addLayout(btn_layout)
        
        # Add the scroll area to the page's main layout
        self.layout().addWidget(scroll)

        # Paths to output files (populated after a run completes)
        self.sps_plot_path = None
        self.gene_ranks_path = None
        self.selected_sites_path = None
    
    def on_enter(self):
        """Update the command display when the page is shown."""
        try:
            # Reset progress indicators
            self.overall_progress_bar.setValue(0)
            self.step_progress_bar.setValue(0)
            self.step_status_label.setText("Ready to run.")
            # Reset run button text
            self.run_btn.setText("Run Analysis")

            # Enable/disable wizard buttons
            if self.wizard():
                self.wizard().button(QWizard.WizardButton.BackButton).setEnabled(True)
                self.wizard().button(QWizard.WizardButton.FinishButton).hide()

            self.run_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)

        except ValueError as e:
            self.cmd_display.setPlainText(f"Error generating command: {str(e)}")
            self.run_btn.setEnabled(False)
    
    def run_analysis(self):
        """Run the ESL-PSC analysis."""
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
        os.chdir(project_root)
        try:
            # Ask for confirmation if output will be overwritten
            output_dir = self.config.output_dir
            base_name = self.config.output_file_base_name
            if os.path.isdir(output_dir):
                # Determine if the current run will overwrite existing files by
                # checking for *exact* filename matches of the outputs that will
                # be produced. 
                expected_output_filenames = [
                    f"{base_name}_gene_ranks.csv",
                    f"{base_name}_species_predictions.csv",
                ]
                existing_conflicts = [
                    fn for fn in expected_output_filenames
                    if os.path.isfile(os.path.join(output_dir, fn))
                ]
                if existing_conflicts:
                    reply = QMessageBox.warning(
                        self,
                        "Overwrite Existing Output?",
                        (
                            "Existing outputs with the same basename, "
                            f"'{base_name}', in the folder '{output_dir}' will "
                            "be overwritten by the outputs from this run.\n\nContinue?"
                        ),
                        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                        QMessageBox.StandardButton.No,
                    )
                    if reply == QMessageBox.StandardButton.No:
                        return

            args = self.config.get_command_args()

            # Update UI for running state
            self.run_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            if self.wizard():
                self.wizard().button(QWizard.WizardButton.BackButton).setEnabled(False)
            
            self.cmd_display.clear()
            self.cmd_display.appendPlainText(f"$ python -m esl_multimatrix.py {self.config.get_command_string()}")
            self.step_status_label.setText("Starting analysis...")
            self.sps_btn.hide()
            self.gene_btn.hide()
            self.sites_btn.hide()
            self.sps_plot_path = None
            self.gene_ranks_path = None
            self.selected_sites_path = None

            # Create and configure worker
            self.worker = ESLWorker(args)
            self.worker.signals.output.connect(self.append_output)
            self.worker.signals.error.connect(self.append_error)
            self.worker.signals.finished.connect(self.analysis_finished)
            self.worker.signals.overall_progress.connect(self.update_overall_progress)
            self.worker.signals.step_progress.connect(self.update_step_progress)
            self.worker.signals.step_status.connect(self.update_step_status)
            
            # Start worker
            self.thread_pool.start(self.worker)
            
        except Exception as e:
            self.append_error(f"Error starting analysis: {str(e)}")
            self.analysis_finished(1)
    
    def stop_analysis(self):
        """Stop the running analysis."""
        if self.worker and self.worker.is_running:
            self.worker.stop()
    
    def append_output(self, text):
        """Append text to the output display."""
        # Using append ensures new lines display correctly in rich text mode
        for line in text.splitlines() or [""]:
            self.cmd_display.appendPlainText(line)

    def append_error(self, text):
        """Append error text to the output display."""
        for line in text.splitlines() or [""]:
            # Prepend 'ERROR' tag to make error lines stand out while keeping plain text formatting
            self.cmd_display.appendPlainText(f"[ERROR] {line}")
    
    def update_overall_progress(self, value):
        """Update the overall progress bar."""
        self.overall_progress_bar.setValue(value)

    def update_step_progress(self, value):
        """Update the step-specific progress bar."""
        self.step_progress_bar.setValue(value)

    def update_step_status(self, text):
        """Update the status label for the current step."""
        self.step_status_label.setText(text)
    
    def show_sps_plot(self):
        """Slot to show the SPS plot if available."""
        if self.sps_plot_path and os.path.exists(self.sps_plot_path):
            SpsPlotDialog.show_dialog(self.sps_plot_path, parent=self)
        else:
            QMessageBox.warning(self, "File Not Found", "The SPS plot file could not be found.")

    def show_gene_ranks(self):
        """Slot to show the gene ranks table if available."""
        if self.gene_ranks_path and os.path.exists(self.gene_ranks_path):
            GeneRanksDialog.show_dialog(
                self.gene_ranks_path,
                self.config,
                self.selected_sites_path,
                parent=self,
            )
        else:
            QMessageBox.warning(self, "File Not Found", "The gene ranks file could not be found.")

    def show_selected_sites(self):
        """Slot to show the selected sites table if available."""
        if self.selected_sites_path and os.path.exists(self.selected_sites_path):
            alignments_dir = getattr(self.config, "alignments_dir", "")
            SelectedSitesDialog.show_dialog(
                self.selected_sites_path,
                self.gene_ranks_path,
                alignments_dir,
                parent=self,
            )
        else:
            QMessageBox.warning(self, "File Not Found", "The selected sites file could not be found.")

    def analysis_finished(self, exit_code):
        """Handle analysis completion."""
        self.run_btn.setEnabled(True)
        self.run_btn.setText("Run New Analysis")
        self.stop_btn.setEnabled(False)

        if self.wizard():
            self.wizard().button(QWizard.WizardButton.BackButton).setEnabled(True)
        
        if exit_code == 0:
            self.overall_progress_bar.setValue(100)
            self.step_progress_bar.setValue(100)
            self.step_status_label.setText("Analysis completed successfully!")
            self.append_output("\n✅ Analysis completed successfully!")

            # Determine output paths and show result buttons
            base = self.config.output_file_base_name
            out_dir = self.config.output_dir

            self.sps_plot_path = None
            if (
                (getattr(self.config, "make_sps_plot", False) or getattr(self.config, "make_sps_kde_plot", False))
                and not getattr(self.config, "no_pred_output", False)
            ):
                fig_name = f"{base}_pred_sps_plot.svg"
                path = os.path.abspath(os.path.join(out_dir, fig_name))
                if os.path.exists(path):
                    self.sps_plot_path = path
                    self.sps_btn.show()

            ranks_path = os.path.abspath(os.path.join(out_dir, f"{base}_gene_ranks.csv"))
            # Only show the button if gene output was requested in this run
            if not getattr(self.config, "no_genes_output", False) and os.path.exists(ranks_path):
                self.gene_ranks_path = ranks_path
                self.gene_btn.show()

            sites_path = os.path.abspath(os.path.join(out_dir, f"{base}_selected_sites.csv"))
            # Only show the button if selected-sites output was generated for this run
            if getattr(self.config, "show_selected_sites", False) and os.path.exists(sites_path):
                self.selected_sites_path = sites_path
                self.sites_btn.show()
        elif exit_code == -1 or (self.worker and self.worker.was_stopped):
            self.step_status_label.setText("Analysis stopped by user.")
            self.append_error("\n🛑 Analysis was stopped.")
        else:
            self.step_status_label.setText("Analysis failed.")
            self.append_error(f"\n❌ Analysis failed with exit code {exit_code}.")
