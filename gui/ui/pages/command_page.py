"""Configuration-review / output page of the ESL-PSC wizard."""
from __future__ import annotations

import os

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QTextEdit, 
    QPushButton, QGroupBox, QFormLayout, QLineEdit, 
    QSizePolicy, QWidget, QScrollArea, QApplication, QFileDialog
)

from .base_page import BaseWizardPage
from .run_page import RunPage

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
        
        # Command display with larger, resizable monospace font
        self.cmd_display = QTextEdit()
        self.cmd_display.setReadOnly(True)
        # Allow selection & copying on macOS and other platforms
        self.cmd_display.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.cmd_display.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse | Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        font = QFont("Courier New")
        font.setPointSize(11)  # Increased font size
        self.cmd_display.setFont(font)
        self.cmd_display.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.cmd_display.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.cmd_display.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Set size policy to allow resizing
        size_policy = self.cmd_display.sizePolicy()
        size_policy.setVerticalPolicy(QSizePolicy.Policy.Expanding)
        size_policy.setHorizontalPolicy(QSizePolicy.Policy.Expanding)
        self.cmd_display.setSizePolicy(size_policy)
        
        # Set minimum size and make it larger
        self.cmd_display.setMinimumSize(600, 300)  # Wider and taller by default
        
        # Copy to clipboard button with fixed size
        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.setMaximumWidth(150)  # Make the button a normal width
        copy_btn.clicked.connect(self.copy_command_to_clipboard)
        
        # Add widgets to command layout with proper spacing
        cmd_layout.addWidget(self.cmd_display)
        cmd_layout.addWidget(copy_btn, alignment=Qt.AlignmentFlag.AlignRight)  # Right-align the button
        cmd_group.setLayout(cmd_layout)
        
        # Add command group to container
        container_layout.addWidget(cmd_group)
        
        # Configuration summary section
        summary_group = QGroupBox("Configuration Summary")
        summary_group.setObjectName("configSummaryGroup")  # Set object name for finding later
        self.summary_group = summary_group 
        summary_layout = QFormLayout()
        summary_group.setLayout(summary_layout)

        # create one monospace QFont that every summary label / value will share
        self.summary_font = QFont("Courier New")
        self.summary_font.setPointSize(10)
        summary_group.setFont(self.summary_font)
        
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
        # Prefix with the full python invocation for display
        full_cmd = f"python -m esl_multimatrix {cmd_str}"
        self.cmd_display.setPlainText(full_cmd)
        
        # Update the configuration summary
        self.update_summary()
    
    def update_summary(self):
        """Update the configuration summary section with organized and formatted information."""
        # Use the reference saved in __init__
        summary_group = self.summary_group
        if summary_group is None:
            return
            
        # Get the layout and clear existing widgets
        layout = summary_group.layout()
        self.clear_layout(layout)
        
        # Set styles for headers
        header_style = """
            font-size: 12px;
            font-weight: bold;
            margin-top: 12px;
            margin-bottom: 6px;
            color: palette(text);
        """
        
        # Add a header for input files
        input_header = QLabel("Input Files")
        input_header.setStyleSheet(header_style)
        layout.addRow(input_header)
        
        # Add input files
        if hasattr(self.config, 'alignments_dir') and self.config.alignments_dir:
            self.add_summary_item(layout, "Alignment Directory:", self.config.alignments_dir)
        if hasattr(self.config, 'species_groups_file') and self.config.species_groups_file:
            self.add_summary_item(layout, "Species Groups:", self.config.species_groups_file)
        if hasattr(self.config, 'species_phenotypes_file') and self.config.species_phenotypes_file:
            self.add_summary_item(layout, "Phenotypes File:", self.config.species_phenotypes_file)
        if hasattr(self.config, 'prediction_alignments_dir') and self.config.prediction_alignments_dir:
            self.add_summary_item(layout, "Prediction Alignments:", self.config.prediction_alignments_dir)
        if hasattr(self.config, 'limited_genes_file') and self.config.limited_genes_file:
            self.add_summary_item(layout, "Limited Genes:", self.config.limited_genes_file)
        if hasattr(self.config, 'response_dir') and self.config.response_dir:
            self.add_summary_item(layout, "Response Directory:", self.config.response_dir)
        
        # Add a header for analysis parameters
        params_header = QLabel("Analysis Parameters")
        params_header.setStyleSheet(header_style)
        layout.addRow(params_header)
        
        # Add hyperparameters
        if hasattr(self.config, 'initial_lambda1') and hasattr(self.config, 'final_lambda1'):
            lambda1_str = f"{self.config.initial_lambda1} → {self.config.final_lambda1}"
            # Show step size for linear grid
            if hasattr(self.config, 'grid_type') and self.config.grid_type == 'linear' and hasattr(self.config, 'num_points') and self.config.num_points:
                lambda1_str += f" (step: {self.config.num_points})"
            self.add_summary_item(layout, "Lambda 1 Range:", lambda1_str)
        
        if hasattr(self.config, 'initial_lambda2') and hasattr(self.config, 'final_lambda2'):
            lambda2_str = f"{self.config.initial_lambda2} → {self.config.final_lambda2}"
            # Show step size for linear grid
            if hasattr(self.config, 'grid_type') and self.config.grid_type == 'linear' and hasattr(self.config, 'num_points') and self.config.num_points:
                lambda2_str += f" (step: {self.config.num_points})"
            self.add_summary_item(layout, "Lambda 2 Range:", lambda2_str)
        
        # Group penalty settings
        if hasattr(self.config, 'group_penalty_type'):
            gp_str = self.config.group_penalty_type.capitalize()
            if (hasattr(self.config, 'initial_gp_value') and 
                hasattr(self.config, 'final_gp_value') and 
                hasattr(self.config, 'gp_step') and 
                self.config.group_penalty_type not in ['median', 'std']):
                gp_str = f"{gp_str}: {self.config.initial_gp_value} → {self.config.final_gp_value} (step: {self.config.gp_step})"
            self.add_summary_item(layout, "Group Penalty:", gp_str)
        
        # Grid type and parameters
        if hasattr(self.config, 'grid_type'):
            grid_info = self.config.grid_type.capitalize()
            if hasattr(self.config, 'num_points'):
                if self.config.grid_type == 'log':
                    grid_info += f" ({self.config.num_points} points)"
                else:  # linear
                    grid_info += f" (step: {self.config.num_points})"
            self.add_summary_item(layout, "Lambda Grid Type:", grid_info)
        
        # Top rank fraction (always show)
        if hasattr(self.config, 'top_rank_frac'):
            self.add_summary_item(layout, "Top Rank Fraction:", str(self.config.top_rank_frac))
        
        # Phenotype names
        if hasattr(self.config, 'pheno_name1') and hasattr(self.config, 'pheno_name2'):
            self.add_summary_item(layout, "Phenotype Comparison:", 
                               f"{self.config.pheno_name1} vs {self.config.pheno_name2}")
        
        # Add a header for output options
        output_header = QLabel("Output Options")
        output_header.setStyleSheet(header_style)
        layout.addRow(output_header)
        
        # Output directory and base name
        if hasattr(self.config, 'output_dir') and self.config.output_dir:
            self.add_summary_item(layout, "Output Directory:", self.config.output_dir)
        
        if hasattr(self.config, 'output_file_base_name') and self.config.output_file_base_name:
            self.add_summary_item(layout, "Output Base Name:", self.config.output_file_base_name)
            
        # Output settings
        output_settings = []
        
        # Keep raw output
        if hasattr(self.config, 'keep_raw_output') and self.config.keep_raw_output:
            output_settings.append("Keep raw output files")
            
        # Show selected sites
        if hasattr(self.config, 'show_selected_sites') and self.config.show_selected_sites:
            output_settings.append("Show selected sites")
            
        # Output type settings – list which output files will actually be generated
        gene_ranks_output = not getattr(self.config, 'no_genes_output', False)
        species_pred_output = not getattr(self.config, 'no_pred_output', False)

        if gene_ranks_output:
            output_settings.append("Gene ranks file")
        if species_pred_output:
            output_settings.append("Species predictions file")

        # If both primary outputs are disabled, make this explicit for the user
        if not gene_ranks_output and not species_pred_output:
            output_settings.append("No primary output files will be saved")
        
        # SPS plot settings – show SPS, KDE, or both depending on selections
        make_sps_plot = getattr(self.config, 'make_sps_plot', False)
        make_kde_plot = getattr(self.config, 'make_sps_kde_plot', False)

        if make_sps_plot and make_kde_plot:
            output_settings.append("Generate SPS and KDE plots")
        elif make_sps_plot:
            output_settings.append("Generate SPS plots")
        elif make_kde_plot:
            output_settings.append("Generate SPS KDE plots")
            
        if output_settings:
            self.add_summary_item(layout, "Output Settings:", "; ".join(output_settings))
            
        # Add a header for advanced options
        advanced_header = QLabel("Advanced Options")
        advanced_header.setStyleSheet(header_style)
        layout.addRow(advanced_header)
        
        # Deletion canceler options
        deletion_text = ""
        cancel_only_partner = getattr(self.config, 'cancel_only_partner', True)
        min_pairs = getattr(self.config, 'min_pairs', 1)
        nix_full_del = getattr(self.config, 'nix_full_deletions', False)

        if cancel_only_partner:
            deletion_text = "Cancel only partner sequences with gaps"
            if min_pairs > 1:
                deletion_text += f" (min pairs: {min_pairs})"
        else:
            deletion_text = "Eliminate any site with one or more gaps"

        if nix_full_del:
            deletion_text += "; Exclude full deletions"

        # Use existing preprocess annotation
        if getattr(self.config, 'use_existing_preprocess', False):
            deletion_text += "; Use existing preprocess"

        self.add_summary_item(layout, "Deletion Handling:", deletion_text)
        if getattr(self.config, 'use_existing_alignments', False) and getattr(self.config, 'canceled_alignments_dir', ''):
            self.add_summary_item(layout, "Existing Canceled Alignments:", self.config.canceled_alignments_dir)
            
        # Null model options
        null_opts = []
        if hasattr(self.config, 'make_null_models') and self.config.make_null_models:
            null_opts.append("Standard null models")
        if hasattr(self.config, 'make_pair_randomized_null_models') and self.config.make_pair_randomized_null_models:
            rand_str = "Pair-randomized null models"
            if hasattr(self.config, 'num_randomized_alignments'):
                rand_str += f" ({self.config.num_randomized_alignments} randomizations)"
            null_opts.append(rand_str)
            
        if null_opts:
            self.add_summary_item(layout, "Null Models:", "; ".join(null_opts))
        
        # Plot options are now consolidated into the output settings above
            
        # Response directory if specified
        if hasattr(self.config, 'response_dir') and self.config.response_dir:
            self.add_summary_item(layout, "Response Directory:", self.config.response_dir)
    
    def add_summary_item(self, layout, label, value):
        """Add a single item to the configuration summary."""
        label_widget = QLabel(label)
        label_widget.setFont(self.summary_font)
        label_widget.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        
        value_widget = QLineEdit(str(value))
        value_widget.setReadOnly(True)
        value_widget.setFont(self.summary_font)
        # No inline stylesheet – rely on global wizard stylesheet so borders adapt to theme.
        value_widget.setMinimumWidth(500)
        
        # Set size policy to allow horizontal expansion
        value_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        
        # Enable horizontal scrolling for long paths
        value_widget.setCursorPosition(0)  # Start scrolled to the beginning
        
        # Create a container widget with horizontal layout to ensure proper expansion
        container = QWidget()
        container_layout = QHBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.addWidget(value_widget)
        
        # Add a stretch to push the field to the left
        container_layout.addStretch()
        
        layout.addRow(label_widget, container)
    
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