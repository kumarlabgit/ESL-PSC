"""Configuration-review / output page of the ESL-PSC wizard."""
from __future__ import annotations

import os

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QFormLayout, QLabel,
    QTextEdit, QPushButton, QLineEdit, QApplication, QFileDialog
)

from .base_page import BaseWizardPage

class CommandPage(BaseWizardPage):
    """Page for reviewing configuration before running the analysis."""
    
    def __init__(self, config):
        """Initialize the output page."""
        super().__init__("Review Configuration")
        self.config = config
        self.setSubTitle("Review your configuration and command before running the analysis.")
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        # Create a container widget for the scroll area
        container = QWidget()
        scroll.setWidget(container)
        
        # Create a layout for the container
        container_layout = QVBoxLayout(container)
        
        # Command preview section
        cmd_group = QGroupBox("ESL-PSC Command")
        cmd_layout = QVBoxLayout()
        
        # Command display with monospace font
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        font = QFont("Monospace")
        font.setStyleHint(QFont.StyleHint.TypeWriter)
        self.cmd_display.setFont(font)
        self.cmd_display.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.cmd_display.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.cmd_display.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Copy to clipboard button
        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(self.copy_command_to_clipboard)
        
        # Add widgets to command layout
        cmd_layout.addWidget(self.cmd_display)
        cmd_layout.addWidget(copy_btn)
        cmd_group.setLayout(cmd_layout)
        
        # Add command group to container
        container_layout.addWidget(cmd_group)
        
        # Configuration summary section
        summary_group = QGroupBox("Configuration Summary")
        summary_layout = QFormLayout()
        summary_group.setLayout(summary_layout)
        
        # Summary fields will be populated in on_enter()
        self.summary_labels = {}
        
        # Add some spacing between fields
        summary_layout.setVerticalSpacing(5)
        
        # Add summary group to container
        container_layout.addWidget(summary_group)
        
        # Add stretch to push everything to the top
        container_layout.addStretch()
        
        # Add the scroll area to the page's layout
        self.layout().addWidget(scroll)
    
    def copy_command_to_clipboard(self):
        """Copy the current command to the system clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.cmd_display.toPlainText())
        
    def browse_output_dir(self):
        """Open a directory dialog to select the output directory."""
        current_dir = self.output_dir_edit.text() or os.getcwd()
        dir_path = QFileDialog.getExistingDirectory(
            self, 
            "Select Output Directory",
            current_dir,
            options=QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks
        )
        
        if dir_path:  # User didn't cancel
            self.config.output_dir = dir_path
            self.output_dir_edit.setText(dir_path)
            # Update the command display to reflect the new output directory
            self.on_enter()
        
    def on_enter(self):
        """Update the command preview and summary when the page is shown."""
        # Use the config's get_command_string method to generate the command
        try:
            cmd_str = self.config.get_command_string()
        except ValueError as e:
            cmd_str = f"Error generating command: {str(e)}"
            
        # Display the command
        self.cmd_display.setPlainText(cmd_str)
        
        # Update the configuration summary
        self.update_summary()
    
    def update_summary(self):
        """Update the configuration summary section."""
        # Get the summary group box
        summary_group = self.findChild(QGroupBox, "Configuration Summary")
        if not summary_group:
            return
            
        # Get the layout and clear existing widgets
        layout = summary_group.layout()
        self.clear_layout(layout)
        
        # Add configuration items
        # Required parameters
        if hasattr(self.config, 'alignments_dir') and self.config.alignments_dir:
            self.add_summary_item(layout, "Alignment Directory:", self.config.alignments_dir)
        if hasattr(self.config, 'species_groups_file') and self.config.species_groups_file:
            self.add_summary_item(layout, "Species Groups File:", self.config.species_groups_file)
        
        # Optional input files
        if hasattr(self.config, 'species_phenotypes_file') and self.config.species_phenotypes_file:
            self.add_summary_item(layout, "Phenotypes File:", self.config.species_phenotypes_file)
        if hasattr(self.config, 'prediction_alignments_dir') and self.config.prediction_alignments_dir:
            self.add_summary_item(layout, "Prediction Alignments Dir:", self.config.prediction_alignments_dir)
        if hasattr(self.config, 'limited_genes_file') and self.config.limited_genes_file:
            self.add_summary_item(layout, "Limited Genes File:", self.config.limited_genes_file)
        
        # Hyperparameters
        if hasattr(self.config, 'initial_lambda1') and hasattr(self.config, 'final_lambda1') and hasattr(self.config, 'lambda_step'):
            lambda1_str = f"{self.config.initial_lambda1} to {self.config.final_lambda1} (step: {self.config.lambda_step})"
            self.add_summary_item(layout, "Lambda 1 Range:", lambda1_str)
        
        if hasattr(self.config, 'initial_lambda2') and hasattr(self.config, 'final_lambda2') and hasattr(self.config, 'lambda2_step'):
            lambda2_str = f"{self.config.initial_lambda2} to {self.config.final_lambda2} (step: {self.config.lambda2_step})"
            self.add_summary_item(layout, "Lambda 2 Range:", lambda2_str)
        
        # Group penalty settings
        if hasattr(self.config, 'group_penalty_type') and hasattr(self.config, 'initial_gp_value') and hasattr(self.config, 'final_gp_value') and hasattr(self.config, 'gp_step'):
            gp_str = f"{self.config.group_penalty_type}: {self.config.initial_gp_value} to {self.config.final_gp_value} (step: {self.config.gp_step})"
            self.add_summary_item(layout, "Group Penalty:", gp_str)
        
        # Logspace and other settings
        if hasattr(self.config, 'use_logspace'):
            logspace_str = f"{'Yes' if self.config.use_logspace else 'No'}"
            if self.config.use_logspace and hasattr(self.config, 'num_log_points'):
                logspace_str += f" ({self.config.num_log_points} points)"
            self.add_summary_item(layout, "Use Logspace:", logspace_str)
        
        # Phenotype names
        if hasattr(self.config, 'pheno_name1') and hasattr(self.config, 'pheno_name2'):
            self.add_summary_item(layout, "Phenotype Names:", f"{self.config.pheno_name1} vs {self.config.pheno_name2}")
        
        # Deletion canceler options
        if hasattr(self.config, 'nix_full_deletions') or hasattr(self.config, 'cancel_only_partner') or hasattr(self.config, 'min_pairs'):
            del_opts = []
            if hasattr(self.config, 'nix_full_deletions') and self.config.nix_full_deletions:
                del_opts.append("Exclude full deletions")
            if hasattr(self.config, 'cancel_only_partner') and not self.config.cancel_only_partner:
                del_opts.append("Cancel all partners")
            if hasattr(self.config, 'min_pairs') and self.config.min_pairs > 1:
                del_opts.append(f"Min pairs: {self.config.min_pairs}")
            
            if del_opts:
                self.add_summary_item(layout, "Deletion Options:", ", ".join(del_opts))
        
        # Output options
        if hasattr(self.config, 'output_dir') and self.config.output_dir:
            self.add_summary_item(layout, "Output Directory:", self.config.output_dir)
        
        if hasattr(self.config, 'output_file_base_name') and self.config.output_file_base_name:
            self.add_summary_item(layout, "Output Base Name:", self.config.output_file_base_name)
        
        # Output toggles
        if hasattr(self.config, 'keep_raw_output') or hasattr(self.config, 'show_selected_sites'):
            toggles = []
            if hasattr(self.config, 'keep_raw_output') and self.config.keep_raw_output:
                toggles.append("Keep raw output")
            if hasattr(self.config, 'show_selected_sites') and self.config.show_selected_sites:
                toggles.append("Show selected sites")
            
            if toggles:
                self.add_summary_item(layout, "Output Options:", ", ".join(toggles))
        
        # Plot options
        if hasattr(self.config, 'make_sps_plot') or hasattr(self.config, 'make_sps_kde_plot'):
            plots = []
            if hasattr(self.config, 'make_sps_plot') and self.config.make_sps_plot:
                plots.append("SPS Plot")
            if hasattr(self.config, 'make_sps_kde_plot') and self.config.make_sps_kde_plot:
                plots.append("SPS KDE Plot")
            
            if plots:
                self.add_summary_item(layout, "Generate Plots:", ", ".join(plots))
        
        # Multi-matrix options
        multi_opts = []
        if hasattr(self.config, 'top_rank_frac') and self.config.top_rank_frac != 0.01:
            multi_opts.append(f"Top rank frac: {self.config.top_rank_frac}")
        if hasattr(self.config, 'response_dir') and self.config.response_dir:
            multi_opts.append(f"Response dir: {self.config.response_dir}")
        
        # Null model options
        if hasattr(self.config, 'make_null_models') and self.config.make_null_models:
            multi_opts.append("Generate null models")
        if hasattr(self.config, 'make_pair_randomized_null_models') and self.config.make_pair_randomized_null_models:
            rand_str = "Generate pair-randomized nulls"
            if hasattr(self.config, 'num_randomized_alignments') and self.config.num_randomized_alignments != 10:
                rand_str += f" ({self.config.num_randomized_alignments} randomizations)"
            multi_opts.append(rand_str)
        
        if multi_opts:
            self.add_summary_item(layout, "Multi-matrix Options:", "; ".join(multi_opts))
    
    def add_summary_item(self, layout, label, value):
        """Add a single item to the configuration summary."""
        label_widget = QLabel(label)
        value_widget = QLineEdit(value)
        value_widget.setReadOnly(True)
        value_widget.setStyleSheet("background: #f8f8f8; border: 1px solid #ddd;")
        layout.addRow(label_widget, value_widget)
    
    def clear_layout(self, layout):
        """Clear all widgets from a layout."""
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    def run_analysis(self):
        """Run the ESL-PSC analysis."""
        # Create a new RunPage instance and show it
        run_page = RunPage(self.config)
        self.wizard().addPage(run_page)
        self.wizard().setCurrentPage(run_page)