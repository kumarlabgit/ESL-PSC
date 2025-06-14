"""Run-analysis page of the ESL-PSC wizard (live output & progress)."""
from __future__ import annotations

from PyQt6.QtCore import Qt, QThreadPool
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QTextEdit, QPushButton,
    QLabel, QProgressBar, QHBoxLayout, QApplication
)

from gui.core.worker import ESLWorker
from .base_page import BaseWizardPage

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
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
        # Command display section
        cmd_group = QGroupBox("ESL-PSC Command")
        cmd_layout = QVBoxLayout()
        
        # Command display with monospace font for better readability
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        self.cmd_display.setFont(QFont("Courier", 10))
        self.cmd_display.setPlaceholderText("The ESL-PSC command will be executed...")
        self.cmd_display.setMinimumHeight(100)
        
        # Copy button
        self.copy_btn = QPushButton("Copy to Clipboard")
        self.copy_btn.clicked.connect(self.copy_command_to_clipboard)
        
        # Add widgets to command layout
        cmd_layout.addWidget(QLabel("Command being executed:"))
        cmd_layout.addWidget(self.cmd_display)
        cmd_layout.addWidget(self.copy_btn, 0, Qt.AlignmentFlag.AlignRight)
        cmd_group.setLayout(cmd_layout)
        container_layout.addWidget(cmd_group)
        
        # Output display
        output_group = QGroupBox("Analysis Output")
        output_layout = QVBoxLayout()
        
        self.output_display = QTextEdit()
        self.output_display.setReadOnly(True)
        self.output_display.setFont(QFont("Courier", 10))
        self.output_display.setPlaceholderText("Analysis output will appear here...")
        output_layout.addWidget(self.output_display)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setTextVisible(True)
        output_layout.addWidget(QLabel("Progress:"))
        output_layout.addWidget(self.progress_bar)
        
        output_group.setLayout(output_layout)
        container_layout.addWidget(output_group, 1)  # Take remaining space
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.go_back)
        
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_analysis)
        
        btn_layout.addWidget(self.back_btn)
        btn_layout.addStretch()
        btn_layout.addWidget(self.run_btn)
        btn_layout.addWidget(self.stop_btn)
        
        container_layout.addLayout(btn_layout)
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)
    
    def copy_command_to_clipboard(self):
        """Copy the current command to the system clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.cmd_display.toPlainText())
    
    def go_back(self):
        """Go back to the previous page."""
        self.wizard().back()
    
    def on_enter(self):
        """Update the command display when the page is shown."""
        try:
            # Get the command string from the config
            cmd_str = self.config.get_command_string()
            
            # Display the command
            self.cmd_display.setPlainText(cmd_str)
            
            # Clear previous output
            self.output_display.clear()
            
            # Reset progress
            self.progress_bar.setValue(0)
            
            # Enable/disable buttons
            self.run_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
            
        except ValueError as e:
            # If there's an error generating the command, show it
            self.cmd_display.setPlainText(f"Error generating command: {str(e)}")
            self.run_btn.setEnabled(False)
            self.stop_btn.setEnabled(False)
    
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
