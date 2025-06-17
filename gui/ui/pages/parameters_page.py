"""Analysis-parameters page of the ESL-PSC wizard."""
from __future__ import annotations

import os

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QWizard
from dataclasses import fields
from gui.core.config import ESLConfig
from PyQt6.QtWidgets import (
    QScrollArea, QWidget, QVBoxLayout, QGroupBox, QFormLayout, QHBoxLayout,
    QLabel, QButtonGroup, QRadioButton, QLineEdit, QDoubleSpinBox, QSpinBox,
    QComboBox, QCheckBox, QPushButton, QStackedWidget, QFrame, QFileDialog,
    QAbstractSpinBox
)

from gui.ui.widgets.file_selectors import FileSelector
from .base_page import BaseWizardPage

class ParametersPage(BaseWizardPage):
    """Page for setting analysis parameters."""
    
    def __init__(self, config, parent=None):
        super().__init__("Analysis Parameters", parent)
        self.config = config
        # Add subtitle to ensure consistent header with divider
        self.setSubTitle("Configure settings and parameters for the analysis and output.")
        self.has_species_pheno = False  # Track if we have a species phenotype file
        
        # Restore Defaults button (will live in the wizard footer)
        self.restore_defaults_btn = QPushButton("Restore Defaults")
        self.restore_defaults_btn.clicked.connect(self.restore_defaults)

        # Store references to widgets that might be accessed after deletion
        self.widgets_initialized = False
        
        # Create scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        
        # Create container widget for scroll area
        container = QWidget()
        self.container_layout = QVBoxLayout(container)
        self.container_layout.setContentsMargins(10, 10, 10, 10)
        self.container_layout.setSpacing(15)
        
        # Add container to scroll area
        scroll.setWidget(container)
        
        # Add scroll area to the page's layout
        layout = self.layout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(scroll)
        
        # ===== Output Options =====
        output_group = QGroupBox("Output Options")
        output_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "Control which output files are generated. Including more outputs will increase "
            "computation time and disk usage."
        )
        explanation.setWordWrap(True)
        output_layout.addWidget(explanation)

        # Output file base name (on its own row, left-aligned)
        output_name_layout = QHBoxLayout()
        output_name_layout.addWidget(QLabel("Output File Base Name:"))
        
        self.output_file_base_name = QLineEdit("esl_psc_results")
        self.output_file_base_name.setMinimumWidth(400)  # Make it wider
        self.output_file_base_name.textChanged.connect(
            lambda t: setattr(self.config, 'output_file_base_name', t)
        )
        # Set the default value in config
        self.config.output_file_base_name = "esl_psc_results"
        output_name_layout.addWidget(self.output_file_base_name)
        output_name_layout.addStretch()  # Push everything to the left
        
        output_layout.addLayout(output_name_layout)
        
        # Reduce spacing before the next section
        output_layout.addSpacing(5)  # Reduced from 10px to 5px

        # Output directory
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(QLabel("Output Directory:"))
        
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setReadOnly(True)
        self.output_dir_edit.setPlaceholderText("Click Browse to select output directory")
        output_dir_layout.addWidget(self.output_dir_edit, 1)  # Allow expanding
        
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse_output_dir)
        output_dir_layout.addWidget(self.browse_btn)
        
        output_layout.addLayout(output_dir_layout)
        
        # Add vertical spacer for better separation
        output_layout.addSpacing(10)  # Add 10px spacing
        
        # Radio buttons for output options
        self.output_options_group = QButtonGroup(self)

        output_files_group = QGroupBox("Output Files to Generate")
        file_layout = QVBoxLayout(output_files_group)
        
        # Option 1: Gene ranks only (fastest)
        genes_only_layout = QHBoxLayout()
        self.genes_only_btn = QRadioButton("Gene ranks only (fastest)")
        self.genes_only_btn.setToolTip(
            "Only generate gene ranks output. This is the fastest option and is "
            "recommended for finding genes that might be related to a convergent phenotype."
        )
        self.genes_only_btn.setChecked(True)  # Set as default
        self.output_options_group.addButton(self.genes_only_btn)
        genes_only_layout.addWidget(self.genes_only_btn)
        genes_only_layout.addStretch()
        file_layout.addLayout(genes_only_layout)
        
        # Option 2: Species predictions only (requires species phenotype file)
        preds_only_layout = QHBoxLayout()
        self.preds_only_btn = QRadioButton("Species predictions only (requires species phenotype file)")
        self.preds_only_btn.setToolTip(
            "Only generate species predictions output. This requires a species phenotype file "
            "to be specified on the input page."
        )
        self.preds_only_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_genes_output', checked)
        )
        self.output_options_group.addButton(self.preds_only_btn)
        preds_only_layout.addWidget(self.preds_only_btn)
        preds_only_layout.addStretch()
        file_layout.addLayout(preds_only_layout)
        
        # Option 3: Both outputs (slowest)
        both_outputs_layout = QHBoxLayout()
        self.both_outputs_btn = QRadioButton("Both outputs (slowest)")
        self.both_outputs_btn.setToolTip(
            "Generate both gene ranks and species predictions outputs. "
            "This is the most comprehensive but slowest option."
        )
        self.both_outputs_btn.setChecked(True)  # Default selection
        self.output_options_group.addButton(self.both_outputs_btn)
        both_outputs_layout.addWidget(self.both_outputs_btn)
        both_outputs_layout.addStretch()
        file_layout.addLayout(both_outputs_layout)

        output_layout.addWidget(output_files_group)
        
        # Connect all radio button signals after they're all created
        self.genes_only_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_pred_output', checked)
        )
        self.genes_only_btn.toggled.connect(self._update_phenotype_names_state)
        self.preds_only_btn.toggled.connect(self._update_phenotype_names_state)
        self.both_outputs_btn.toggled.connect(self._update_phenotype_names_state)
        
        # Add more spacing after the output options
        output_layout.addSpacing(15)  # Increased spacing
        
        # Phenotype Names
        # Add label above the fields
        pheno_label = QLabel("Phenotype Names (for output files):")
        output_layout.addWidget(pheno_label)
        
        # Create a container widget for the fields
        pheno_container = QWidget()
        pheno_layout = QHBoxLayout(pheno_container)
        pheno_layout.setContentsMargins(0, 0, 0, 0)
        pheno_layout.setSpacing(10)
        
        # Create fields with labels
        pos_container = QWidget()
        pos_layout = QHBoxLayout(pos_container)
        pos_layout.setContentsMargins(0, 0, 0, 0)
        pos_layout.addWidget(QLabel("Positive:"))
        
        self.pheno_name1 = QLineEdit("Convergent")
        self.pheno_name1.setPlaceholderText("Positive phenotype name")
        self.pheno_name1.textChanged.connect(
            lambda t: setattr(self.config, 'pheno_name1', t)
        )
        pos_layout.addWidget(self.pheno_name1)
        
        neg_container = QWidget()
        neg_layout = QHBoxLayout(neg_container)
        neg_layout.setContentsMargins(10, 0, 0, 0)
        neg_layout.addWidget(QLabel("Negative:"))
        
        self.pheno_name2 = QLineEdit("Control")
        self.pheno_name2.setPlaceholderText("Negative phenotype name")
        self.pheno_name2.textChanged.connect(
            lambda t: setattr(self.config, 'pheno_name2', t)
        )
        neg_layout.addWidget(self.pheno_name2)
        
        # Set default values in config
        self.config.pheno_name1 = "Convergent"
        self.config.pheno_name2 = "Control"
        
        # Add fields to layout
        pheno_layout.addWidget(pos_container)
        pheno_layout.addWidget(neg_container)
        pheno_layout.addStretch()
        
        # Add to main layout
        output_layout.addWidget(pheno_container)
        
        # Store references for enabling/disabling
        self.pheno_label = pheno_label
        self.pheno_container = pheno_container
        
        # Initialize the state of phenotype names based on output option
        self._update_phenotype_names_state()
        
        # Add standard spacing after phenotype names
        output_layout.addSpacing(8)
        
        # Add vertical spacer for better separation before output directory
        output_layout.addSpacing(5)  # Reduced spacing before output directory
        
        # SPS plot options (always visible but may be disabled)
        self.sps_plot_group = QGroupBox("Species Prediction Score (SPS) Plots")
        sps_plot_layout = QVBoxLayout()
        sps_plot_layout.setSpacing(5)  # Reduce spacing between radio buttons
        
        # Create a button group for mutually exclusive plot options
        self.plot_options_group = QButtonGroup(self)
        
        # None option (default)
        self.no_sps_plot = QRadioButton("None")
        self.no_sps_plot.setToolTip("Do not generate any SPS plots.")
        self.no_sps_plot.setChecked(True)  # Default selection
        self.no_sps_plot.toggled.connect(
            lambda checked: setattr(self.config, 'no_sps_plot', checked)
        )
        self.plot_options_group.addButton(self.no_sps_plot)
        sps_plot_layout.addWidget(self.no_sps_plot)
        
        # Make SPS plot
        self.make_sps_plot = QRadioButton("Generate SPS Violin plots")
        self.make_sps_plot.setToolTip(
            "Create violin plots showing SPS density for each true phenotype."
        )
        self.make_sps_plot.toggled.connect(
            lambda checked: setattr(self.config, 'make_sps_plot', checked)
        )
        self.plot_options_group.addButton(self.make_sps_plot)
        sps_plot_layout.addWidget(self.make_sps_plot)
        
        # Make SPS KDE plot
        self.make_sps_kde_plot = QRadioButton("Generate SPS KDE plots")
        self.make_sps_kde_plot.setToolTip(
            "Create Kernel Density Estimate (KDE) plots showing SPS density for each true phenotype."
        )
        self.make_sps_kde_plot.toggled.connect(
            lambda checked: setattr(self.config, 'make_sps_kde_plot', checked)
        )
        self.plot_options_group.addButton(self.make_sps_kde_plot)
        sps_plot_layout.addWidget(self.make_sps_kde_plot)
        
        # Set initial config values
        self.config.no_sps_plot = True
        self.config.make_sps_plot = False
        self.config.make_sps_kde_plot = False
        
        self.sps_plot_group.setLayout(sps_plot_layout)
        output_layout.addWidget(self.sps_plot_group)
        
        # Connect output type changes to enable/disable SPS plot options
        def update_sps_plot_state():
            enable_sps = not self.genes_only_btn.isChecked()
            # Set the entire group box and its children enabled/disabled
            self.sps_plot_group.setEnabled(enable_sps)
            # Visually indicate the disabled state
            self.sps_plot_group.setStyleSheet(
                "QGroupBox:disabled { color: gray; }"
                "QCheckBox:disabled { color: gray; }"
            )
        
        self.genes_only_btn.toggled.connect(update_sps_plot_state)
        self.preds_only_btn.toggled.connect(update_sps_plot_state)
        self.both_outputs_btn.toggled.connect(update_sps_plot_state)
        update_sps_plot_state()  # Initial update

        # Keep raw output
        self.keep_raw_output_chk = QCheckBox("Keep raw model weights files")
        self.keep_raw_output_chk.setToolTip("If checked, keep model weights files generated during analysis.")
        self.keep_raw_output_chk.stateChanged.connect(
            lambda s: setattr(self.config, 'keep_raw_output', s == 2)  # 2 is Qt.Checked
        )
        # Reflect existing state if re-entering the page
        self.keep_raw_output_chk.setChecked(getattr(self.config, 'keep_raw_output', False))
        output_layout.addWidget(self.keep_raw_output_chk)
        
        # Show selected sites
        self.show_selected_sites = QCheckBox("Show selected sites in output (currently slow)")
        self.show_selected_sites.setToolTip(
            "If checked, include a dictionary of all selected sites with their highest model score for every gene."
        )
        self.show_selected_sites.stateChanged.connect(
            lambda s: setattr(self.config, 'show_selected_sites', s == 2)
        )
        output_layout.addWidget(self.show_selected_sites)
        
        # Add vertical spacer for better separation
        output_layout.addSpacing(10)  # Add 10px spacing
        
        output_group.setLayout(output_layout)
        # Add output group to container
        self.container_layout.addWidget(output_group)
        
        # ===== Deletion Canceler Options =====
        del_cancel_group = QGroupBox("Deletion Canceler Options")
        del_cancel_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns with gaps "
            "that could unbalance the analysis."
        )
        explanation.setWordWrap(True)
        del_cancel_layout.addWidget(explanation)
        
        # Nix fully canceled genes
        self.nix_full_deletions = QCheckBox("Exclude fully canceled genes")
        self.nix_full_deletions.setToolTip(
            "If checked, genes that are fully canceled in enough species will be excluded from analysis. "
            "This is equivalent to the --nix_full_deletions command line option."
        )
        self.nix_full_deletions.stateChanged.connect(
            lambda s: setattr(self.config, 'nix_full_deletions', s == 2)  # 2 is Qt.Checked
        )
        del_cancel_layout.addWidget(self.nix_full_deletions)
        
        # Cancel only partner
        self.cancel_only_partner = QCheckBox("Only cancel partner deletions")
        self.cancel_only_partner.setToolTip(
            "If checked, only cancel out the partner of a sequence with a deletion at a site. "
            "If unchecked, sites with any gaps will be fully eliminated."
        )
        self.cancel_only_partner.setChecked(True)  # Set checked by default
        self.cancel_only_partner.stateChanged.connect(
            self._update_deletion_canceler_state
        )
        del_cancel_layout.addWidget(self.cancel_only_partner)
        
        # Minimum aligned pairs
        min_pairs_layout = QHBoxLayout()
        min_pairs_layout.addWidget(QLabel("Minimum aligned pairs:"))
        self.min_pairs = QSpinBox()
        self.min_pairs.setRange(2, 100)  # Minimum is 2 as per requirements
        self.min_pairs.setValue(2)
        self.min_pairs.setEnabled(True)  # Enabled by default since checkbox is checked
        self.min_pairs.setToolTip(
            "Minimum number of ungapped aligned pairs required to consider a site. "
            "Sites with fewer ungapped aligned pairs will be excluded from analysis. "

        )
        self.min_pairs.valueChanged.connect(
            lambda v: setattr(self.config, 'min_pairs', v)
        )
        
        # Set initial config values
        self.config.cancel_only_partner = True
        self.config.min_pairs = 2
        min_pairs_layout.addWidget(self.min_pairs)
        min_pairs_layout.addStretch()
        del_cancel_layout.addLayout(min_pairs_layout)
        
        del_cancel_group.setLayout(del_cancel_layout)
        self.container_layout.addWidget(del_cancel_group)

        # ===== Hyperparameters Section =====
        hyper_group = QGroupBox("Hyperparameters")
        hyper_layout = QFormLayout()
        # Ensure the form rows themselves stay left-aligned instead of centering in the
        # available space.  Using AlignLeft keeps the widgets flush with the left
        # edge without forcing them to expand across the full width.
        hyper_layout.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        
        # Grid Type Selection
        grid_type_group = QHBoxLayout()
        self.grid_type_btns = QButtonGroup()
        
        self.logspace_btn = QRadioButton("Logarithmic Grid (Recommended)")
        self.linear_btn = QRadioButton("Linear Grid")
        
        self.grid_type_btns.addButton(self.logspace_btn, 0)
        self.grid_type_btns.addButton(self.linear_btn, 1)
        
        # Set default grid type
        self.config.grid_type = 'log'  # Default to log grid
        self.logspace_btn.setChecked(True)
        
        # Connect signals to update config
        self.logspace_btn.toggled.connect(
            lambda checked: setattr(self.config, 'grid_type', 'log' if checked else 'linear')
        )
        
        grid_type_group.addWidget(self.logspace_btn)
        grid_type_group.addWidget(self.linear_btn)
        grid_type_group.addStretch()
        
        hyper_layout.addRow("Lambda Grid Type:", grid_type_group)
        
        # Lambda 1 (Position Sparsity)
        lambda1_group = QVBoxLayout()
        lambda1_range_group = QHBoxLayout()
        
        self.initial_lambda1 = QDoubleSpinBox()
        self.initial_lambda1.setRange(0.001, 0.999)
        self.initial_lambda1.setDecimals(3)
        self.initial_lambda1.setSingleStep(0.01)
        self.initial_lambda1.setValue(0.01)
        self.initial_lambda1.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_lambda1', v)
        )
        
        self.final_lambda1 = QDoubleSpinBox()
        self.final_lambda1.setRange(0.001, 0.999)
        self.final_lambda1.setDecimals(3)
        self.final_lambda1.setSingleStep(0.01)
        self.final_lambda1.setValue(0.99)
        self.final_lambda1.valueChanged.connect(
            lambda v: setattr(self.config, 'final_lambda1', v)
        )
        
        # Add range controls to the range group
        lambda1_range_group.addWidget(QLabel("From:"))
        lambda1_range_group.addWidget(self.initial_lambda1)
        lambda1_range_group.addWidget(QLabel("To:"))
        lambda1_range_group.addWidget(self.final_lambda1)
        lambda1_range_group.addStretch()
        
        # Create widgets for step size (linear grid)
        step_widget = QWidget()
        step_layout = QHBoxLayout(step_widget)
        step_layout.setContentsMargins(0, 0, 0, 0)
        
        self.lambda_step = QDoubleSpinBox()
        self.lambda_step.setRange(0.01, 1.0)
        self.lambda_step.setSingleStep(0.01)
        self.lambda_step.setValue(0.1)  # Default step size for linear grid
        self.lambda_step.setDecimals(2)  # Show 2 decimal places
        
        step_layout.addWidget(QLabel("Step:"))
        step_layout.addWidget(self.lambda_step)
        step_layout.addStretch()
        
        # Create widgets for log points (log grid)
        log_widget = QWidget()
        log_layout = QHBoxLayout(log_widget)
        log_layout.setContentsMargins(0, 0, 0, 0)
        
        self.num_log_points = QSpinBox()
        self.num_log_points.setRange(4, 1000)  # Minimum 4 points for log grid
        self.num_log_points.setValue(20)  # Default number of points
        self.num_log_points.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.UpDownArrows)
        
        log_layout.addWidget(QLabel("Points:"))
        log_layout.addWidget(self.num_log_points)
        log_layout.addStretch()
        
        # Set initial config values
        self.config.num_points = 20  # Default number of points for log grid
        
        # Create stacked widgets for lambda1 and lambda2
        self.lambda1_stack = QStackedWidget()
        self.lambda1_stack.addWidget(step_widget)
        self.lambda1_stack.addWidget(log_widget)
        
        # Add both groups to the main lambda1 group
        lambda1_group.addLayout(lambda1_range_group)
        lambda1_group.addWidget(self.lambda1_stack)
        
        # Lambda 2 (Group Sparsity)
        lambda2_group = QVBoxLayout()
        lambda2_range_group = QHBoxLayout()
        
        self.initial_lambda2 = QDoubleSpinBox()
        self.initial_lambda2.setRange(0.001, 0.999)
        self.initial_lambda2.setDecimals(3)
        self.initial_lambda2.setSingleStep(0.01)
        self.initial_lambda2.setValue(0.01)
        self.initial_lambda2.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_lambda2', v)
        )
        
        self.final_lambda2 = QDoubleSpinBox()
        self.final_lambda2.setRange(0.001, 0.999)
        self.final_lambda2.setDecimals(3)
        self.final_lambda2.setSingleStep(0.01)
        self.final_lambda2.setValue(0.99)
        self.final_lambda2.valueChanged.connect(
            lambda v: setattr(self.config, 'final_lambda2', v)
        )
        
        # Add range controls to the range group
        lambda2_range_group.addWidget(QLabel("From:"))
        lambda2_range_group.addWidget(self.initial_lambda2)
        lambda2_range_group.addWidget(QLabel("To:"))
        lambda2_range_group.addWidget(self.final_lambda2)
        lambda2_range_group.addStretch()
        
        # Create a new stacked widget for lambda2
        self.lambda2_stack = QStackedWidget()
        self.lambda2_stack.addWidget(step_widget)
        self.lambda2_stack.addWidget(log_widget)
        
        # Add both groups to the main lambda2 group
        lambda2_group.addLayout(lambda2_range_group)
        lambda2_group.addWidget(self.lambda2_stack)
        
        # Initially set the correct view based on grid type
        self._update_grid_type_view()
        
        # Connect grid type change to update the stacks and config
        self.logspace_btn.toggled.connect(self._update_grid_type_view)
        self.linear_btn.toggled.connect(self._update_grid_type_view)
        
        # Connect spinbox changes to update config
        self.num_log_points.valueChanged.connect(
            lambda v: setattr(self.config, 'num_points', int(v)) if self.logspace_btn.isChecked() else None
        )
        self.lambda_step.valueChanged.connect(
            lambda v: setattr(self.config, 'num_points', float(v)) if self.linear_btn.isChecked() else None
        )
        
        # Set initial values
        self.num_log_points.setValue(20)  # Default for log grid
        self.lambda_step.setValue(0.1)  # Default for linear grid
        
        # Add to the main layout
        hyper_layout.addRow("Lambda 1 (Position Sparsity):", lambda1_group)
        hyper_layout.addRow("Lambda 2 (Group Sparsity):", lambda2_group)
        
        # Group Penalty Settings
        penalty_group = QVBoxLayout()
        
        # Penalty type selection
        type_group = QHBoxLayout()
        self.group_penalty_type = QComboBox()
        self.group_penalty_type.addItems(["median (Recommended)", "std", "sqrt", "linear"])
        self.group_penalty_type.setCurrentText("median (Recommended)")
        self.group_penalty_type.currentTextChanged.connect(
            lambda t: self._update_penalty_type(t.split(' ')[0])  # Remove " (Recommended)" from text
        )
        
        # Set initial group penalty type in config
        self.config.group_penalty_type = "median"
        
        type_group.addWidget(QLabel("Type:"))
        type_group.addWidget(self.group_penalty_type)
        type_group.addStretch()
        
        penalty_group.addLayout(type_group)
        
        # Range settings (will be hidden for median/standard)
        self.penalty_range_group = QWidget()
        range_layout = QHBoxLayout(self.penalty_range_group)
        
        # Set up initial group penalty values
        self.initial_gp_value = QDoubleSpinBox()
        self.initial_gp_value.setRange(0.0, 100.0)
        self.initial_gp_value.setSingleStep(0.5)
        self.initial_gp_value.setValue(1.0)  # Default value for standard
        self.initial_gp_value.valueChanged.connect(
            lambda v: setattr(self.config, 'initial_gp_value', v)
        )
        
        self.final_gp_value = QDoubleSpinBox()
        self.final_gp_value.setRange(0.0, 100.0)
        self.final_gp_value.setSingleStep(0.5)
        self.final_gp_value.setValue(1.0)  # Default value for standard
        self.final_gp_value.valueChanged.connect(
            lambda v: setattr(self.config, 'final_gp_value', v)
        )
        
        self.gp_step = QDoubleSpinBox()
        self.gp_step.setRange(0.1, 10.0)
        self.gp_step.setSingleStep(0.1)
        self.gp_step.setValue(0.1)
        self.gp_step.valueChanged.connect(
            lambda v: setattr(self.config, 'gp_step', v)
        )
        
        range_layout.addWidget(QLabel("From:"))
        range_layout.addWidget(self.initial_gp_value)
        range_layout.addWidget(QLabel("To:"))
        range_layout.addWidget(self.final_gp_value)
        range_layout.addWidget(QLabel("Step:"))
        range_layout.addWidget(self.gp_step)
        
        # Set initial config values before adding to layout
        self.config.initial_gp_value = 1.0
        self.config.final_gp_value = 1.0
        self.config.gp_step = 0.1
        
        # Add the penalty group to the layout
        hyper_layout.addRow("Group Penalty Settings:", penalty_group)
        
        # Add the range group to the penalty group
        penalty_group.addWidget(self.penalty_range_group)
        
        # Initialize the penalty type UI state after all widgets are in place
        self._update_penalty_type("median")
        
        # Add Top rank fraction to Hyperparameters section
        self.top_rank_frac = QDoubleSpinBox()
        self.top_rank_frac.setRange(0.0001, 1.0)
        self.top_rank_frac.setSingleStep(0.01)
        # Ensure the widget doesn't expand excessively on some platforms (e.g., Linux)
        self.top_rank_frac.setDecimals(4)
        self.top_rank_frac.setMaximumWidth(120)
        self.top_rank_frac.setValue(0.01)
        self.top_rank_frac.setToolTip(
            "Fraction of genes to consider top-ranked for the purpose of ranking genes in multimatrix runs. "
        )
        self.top_rank_frac.valueChanged.connect(
            lambda v: setattr(self.config, 'top_rank_frac', v)
        )
        hyper_layout.addRow("Top rank fraction:", self.top_rank_frac)
        
        hyper_group.setLayout(hyper_layout)
        # Add hyper group to container
        self.container_layout.addWidget(hyper_group)
        
        # ===== Null Models Section =====
        null_models_group = QGroupBox("Null Models")
        null_models_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "These methods can be used to generate null models in which the convergent signal is scrambled. "
            "Note that you will get *only* null models, so you should do a separate run first to generate non-null models."
        )
        explanation.setWordWrap(True)
        null_models_layout.addWidget(explanation)
        
        # Radio buttons for null model selection
        self.null_models_group = QButtonGroup(self)
        
        # No null models (fastest)
        no_null_layout = QHBoxLayout()
        self.no_null_btn = QRadioButton("No null models")
        self.no_null_btn.setToolTip(
            "Do not generate null models. Real models will be generated. "
        )
        self.no_null_btn.toggled.connect(
            lambda checked: setattr(self.config, 'no_null_models', checked)
        )
        self.null_models_group.addButton(self.no_null_btn)
        no_null_layout.addWidget(self.no_null_btn)
        no_null_layout.addStretch()
        null_models_layout.addLayout(no_null_layout)
        
        # Response-flipped null models
        response_flip_layout = QHBoxLayout()
        self.response_flip_btn = QRadioButton("Response-flipped null models")
        self.response_flip_btn.setToolTip(
            "Generate null models by randomly flipping response values. "
        )
        self.response_flip_btn.toggled.connect(
            lambda checked: setattr(self.config, 'make_null_models', checked)
        )
        self.null_models_group.addButton(self.response_flip_btn)
        response_flip_layout.addWidget(self.response_flip_btn)
        response_flip_layout.addStretch()
        null_models_layout.addLayout(response_flip_layout)
        
        # Pair-randomized null models
        pair_randomized_layout = QHBoxLayout()
        self.pair_randomized_btn = QRadioButton("Pair-randomized null models")
        self.pair_randomized_btn.setToolTip(
            "Generate null models by randomly flipping or not flipping the AAs between each pair at each site."
        )
        self.pair_randomized_btn.toggled.connect(
            lambda checked: (
                setattr(self.config, 'make_pair_randomized_null_models', checked) or True,
                update_num_rand_visibility()
            )
        )
        self.null_models_group.addButton(self.pair_randomized_btn)
        pair_randomized_layout.addWidget(self.pair_randomized_btn)
        pair_randomized_layout.addStretch()
        null_models_layout.addLayout(pair_randomized_layout)
        
        # Number of randomizations
        num_rand_layout = QHBoxLayout()
        self.num_rand_label = QLabel("Number of randomizations:")
        self.num_rand = QSpinBox()
        self.num_rand.setRange(1, 1000)
        self.num_rand.setValue(10)
        self.num_rand.setToolTip(
            "Number of pair-randomized alignments to generate for the null model."
        )
        self.num_rand.valueChanged.connect(
            lambda v: setattr(self.config, 'num_randomized_alignments', v)
        )
        num_rand_layout.addWidget(self.num_rand_label)
        num_rand_layout.addWidget(self.num_rand)
        num_rand_layout.addStretch()
        
        # Function to update visibility and enabled state of num_rand based on selection
        def update_num_rand_visibility():
            # Only show and enable if pair-randomized null models are selected
            enabled = self.pair_randomized_btn.isChecked()
            self.num_rand_label.setVisible(enabled)
            self.num_rand.setVisible(enabled)
            self.num_rand.setEnabled(enabled)
        
        # Connect signals
        self.no_null_btn.toggled.connect(update_num_rand_visibility)
        self.response_flip_btn.toggled.connect(update_num_rand_visibility)
        self.pair_randomized_btn.toggled.connect(update_num_rand_visibility)
        
        # Initial update
        update_num_rand_visibility()
        
        null_models_layout.addLayout(num_rand_layout)
        
        # Set default selection
        self.no_null_btn.setChecked(True)
        
        null_models_group.setLayout(null_models_layout)
        self.container_layout.addWidget(null_models_group)
        # --- Restore Defaults button (placed at bottom of Parameters page) ---
        self.container_layout.addSpacing(12)
        # Create a horizontal layout for the button to prevent it from stretching
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.restore_defaults_btn)
        button_layout.addStretch()  # Push button to the left
        self.container_layout.addLayout(button_layout)
        
        # Mark widgets as initialized
        self.widgets_initialized = True
        
    def restore_defaults(self):
        """Reset all parameters to their default ESLConfig values and update the UI."""
        defaults = ESLConfig()
        # Do NOT reset input file paths or other pages' settings
        skip_fields = {
            "alignments_dir", "species_groups_file", "species_phenotypes_file",
            "prediction_alignments_dir", "limited_genes_file", "response_dir"
        }
        for f in fields(ESLConfig):
            if f.name in skip_fields:
                continue
            setattr(self.config, f.name, getattr(defaults, f.name))

        # ─── Update widgets ────────────────────────────────────────────────────
        # Grid type & lambda grids
        self.logspace_btn.setChecked(self.config.grid_type == 'log')
        self.linear_btn.setChecked(self.config.grid_type == 'linear')
        self.initial_lambda1.setValue(self.config.initial_lambda1)
        self.final_lambda1.setValue(self.config.final_lambda1)
        self.initial_lambda2.setValue(self.config.initial_lambda2)
        self.final_lambda2.setValue(self.config.final_lambda2)
        self.lambda_step.setValue(self.config.lambda_step)
        self.num_log_points.setValue(self.config.num_points)

        # Group penalty
        self.group_penalty_type.setCurrentText("median (Recommended)")
        self.initial_gp_value.setValue(self.config.initial_gp_value)
        self.final_gp_value.setValue(self.config.final_gp_value)
        self.gp_step.setValue(self.config.gp_step)

        # Top-rank frac
        self.top_rank_frac.setValue(self.config.top_rank_frac)

        # Output options
        self.genes_only_btn.setChecked(False)
        self.preds_only_btn.setChecked(False)
        self.both_outputs_btn.setChecked(True)

        self.output_file_base_name.setText(self.config.output_file_base_name)
        self.pheno_name1.setText(self.config.pheno_name1)
        self.pheno_name2.setText(self.config.pheno_name2)
        self.output_dir_edit.setText(self.config.output_dir)

        self.keep_raw_output_chk.setChecked(self.config.keep_raw_output)
        self.show_selected_sites.setChecked(self.config.show_selected_sites)

        # SPS plot radios
        self.no_sps_plot.setChecked(True)

        # Deletion-canceler
        self.nix_full_deletions.setChecked(self.config.nix_full_deletions)
        self.cancel_only_partner.setChecked(self.config.cancel_only_partner)
        self.min_pairs.setValue(self.config.min_pairs)

        # Null models
        self.no_null_btn.setChecked(True)

        # Trigger dependent UI updates
        self._update_grid_type_view()
        self._update_penalty_type(self.config.group_penalty_type)
        self._update_phenotype_names_state()
        self.update_output_options_state()

    # ──────────────────────────────────────────────────────────────────────────
    # Public helper: sync widgets when configuration is loaded
    # ──────────────────────────────────────────────────────────────────────────
    def update_ui_from_config(self):
        """Synchronize all parameter widgets with values from self.config."""
        if not getattr(self, 'widgets_initialized', False):
            return  # widgets not ready yet
        cfg = self.config
        # Grid type & lambda values
        self.logspace_btn.setChecked(cfg.grid_type == 'log')
        self.linear_btn.setChecked(cfg.grid_type == 'linear')
        self.initial_lambda1.setValue(cfg.initial_lambda1)
        self.final_lambda1.setValue(cfg.final_lambda1)
        self.initial_lambda2.setValue(cfg.initial_lambda2)
        self.final_lambda2.setValue(cfg.final_lambda2)
        if cfg.grid_type == 'log':
            self.num_log_points.setValue(int(cfg.num_points))
        else:
            self.lambda_step.setValue(float(cfg.num_points))
        # Group penalty
        mapping = {
            'median': 'median (Recommended)',
            'standard': 'standard',
            'sqrt': 'sqrt',
            'linear': 'linear'
        }
        self.group_penalty_type.setCurrentText(mapping.get(cfg.group_penalty_type, cfg.group_penalty_type))
        self.initial_gp_value.setValue(cfg.initial_gp_value)
        self.final_gp_value.setValue(cfg.final_gp_value)
        self.gp_step.setValue(cfg.gp_step)
        # Top rank frac
        self.top_rank_frac.setValue(cfg.top_rank_frac)
        # Output basics
        self.output_file_base_name.setText(cfg.output_file_base_name)
        self.pheno_name1.setText(cfg.pheno_name1)
        self.pheno_name2.setText(cfg.pheno_name2)
        self.output_dir_edit.setText(cfg.output_dir)
        # Toggles
        self.keep_raw_output_chk.setChecked(cfg.keep_raw_output)
        self.show_selected_sites.setChecked(cfg.show_selected_sites)
        # Plot radios if exist
        if hasattr(self, 'no_sps_plot'):
            self.no_sps_plot.setChecked(cfg.no_sps_plot)
        if hasattr(self, 'make_sps_plot'):
            self.make_sps_plot.setChecked(cfg.make_sps_plot)
        if hasattr(self, 'make_sps_kde_plot'):
            self.make_sps_kde_plot.setChecked(cfg.make_sps_kde_plot)
        # Deletion canceler
        self.nix_full_deletions.setChecked(cfg.nix_full_deletions)
        self.cancel_only_partner.setChecked(cfg.cancel_only_partner)
        self.min_pairs.setValue(cfg.min_pairs)
        # Update Output Files radios from config (block signals so we don’t write back immediately)
        for btn in (self.genes_only_btn, self.preds_only_btn, self.both_outputs_btn):
            btn.blockSignals(True)

        if cfg.no_genes_output and not cfg.no_pred_output:
            # only species predictions
            self.preds_only_btn.setChecked(True)
        elif cfg.no_pred_output and not cfg.no_genes_output:
            # only gene ranks
            self.genes_only_btn.setChecked(True)
        else:
            # both outputs
            self.both_outputs_btn.setChecked(True)

        for btn in (self.genes_only_btn, self.preds_only_btn, self.both_outputs_btn):
            btn.blockSignals(False)
        # Null models
        self.no_null_btn.setChecked(not (cfg.make_null_models or cfg.make_pair_randomized_null_models))
        self.response_flip_btn.setChecked(cfg.make_null_models)
        self.pair_randomized_btn.setChecked(cfg.make_pair_randomized_null_models)
        self.num_rand.setValue(cfg.num_randomized_alignments)
        # Trigger dependent updates
        self._update_grid_type_view()
        self._update_penalty_type(cfg.group_penalty_type)
        self._update_phenotype_names_state()
        self.update_output_options_state()

    # Qt calls this each time the page becomes current
    def initializePage(self):
        super().initializePage()
        wiz: QWizard = self.wizard()
        if wiz is None:
            return

        # Update UI that depends on input page
        if hasattr(self.wizard(), 'input_page') and hasattr(self.wizard().input_page, 'species_phenotypes'):
            self.has_species_pheno = bool(self.wizard().input_page.species_phenotypes.get_path())
        self.update_output_options_state()
        return

    def _on_wizard_page_changed(self, page_id):
        wiz: QWizard = self.wizard()
        if wiz.currentPage() is self:
            return  # Still on this page
        # Disconnect this slot to avoid repeated calls
        try:
            wiz.currentIdChanged.disconnect(self._on_wizard_page_changed)
        except TypeError:
            pass
        # Hide/disable button and restore old layout
        if wiz.button(QWizard.WizardButton.CustomButton1):
            wiz.button(QWizard.WizardButton.CustomButton1).hide()
        wiz.setButton(QWizard.WizardButton.CustomButton1, None)
        wiz.setOption(QWizard.WizardOption.HaveCustomButton1, False)
        # Restore default button layout so Restore button slot disappears
        if hasattr(self, "_default_button_layout"):
            wiz.setButtonLayout(self._default_button_layout)
        if getattr(self, "_old_button_layout", None):
             wiz.setButtonLayout(self._old_button_layout)
             del self._old_button_layout

    # Also handle Qt Back navigation (cleanupPage called only on Back)
    def cleanupPage(self):
        """Called by QWizard when leaving the page. No per-page cleanup needed."""
        super().cleanupPage()
    
    def browse_output_dir(self):
        """Open a dialog to select the output directory."""
        try:
            # Get the current output directory or use the default
            current_dir = getattr(self.config, 'output_dir', os.getcwd())
            
            # Open directory dialog
            dir_path = QFileDialog.getExistingDirectory(
                self,
                "Select Output Directory",
                current_dir,
                QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks
            )
            
            if dir_path:  # User didn't cancel
                self.output_dir_edit.setText(dir_path)
                self.config.output_dir = dir_path
                
        except Exception as e:
            print(f"Error browsing for output directory: {e}")
            
    def update_output_options_state(self):
        """Update radio-button availability after the phenotype file changes or a config is loaded."""
        # Always re-evaluate whether a phenotype file is present
        try:
            wiz = self.wizard()
            if wiz and hasattr(wiz, 'input_page'):
                self.has_species_pheno = bool(wiz.input_page.species_phenotypes.get_path())
            else:
                # Fallback to the value already stored in the config
                self.has_species_pheno = bool(getattr(self.config, 'species_phenotypes_file', ''))
        except Exception:
            self.has_species_pheno = bool(getattr(self.config, 'species_phenotypes_file', ''))

        if not getattr(self, 'widgets_initialized', False):
            return
            
        try:
            # Check if widgets still exist
            if not hasattr(self, 'preds_only_btn') or not self.preds_only_btn:
                return
                
            # Update UI state based on our flag
            try:
                # Only update if state would change
                if self.preds_only_btn.isEnabled() != self.has_species_pheno:
                    self.preds_only_btn.setEnabled(self.has_species_pheno)
                    self.both_outputs_btn.setEnabled(self.has_species_pheno)
                
                # Update tooltips
                if self.has_species_pheno:
                    self.preds_only_btn.setToolTip(
                        "Only generate species predictions output. This requires a species phenotype file "
                        "to be specified on the input page."
                    )
                    self.both_outputs_btn.setToolTip(
                        "Generate both gene ranks and species predictions outputs. "
                        "This is the most comprehensive but slowest option."
                    )
                else:
                    self.preds_only_btn.setToolTip(
                        "<font color='red'>Requires species phenotype file to be specified on the input page.</font>"
                    )
                    self.both_outputs_btn.setToolTip(
                        "<font color='red'>Requires species phenotype file to be specified on the input page.</font>"
                    )
                
                # If species predictions were selected but no phenotype file, switch to gene ranks only
                if not self.has_species_pheno and (self.preds_only_btn.isChecked() or self.both_outputs_btn.isChecked()):
                    self.genes_only_btn.setChecked(True)
                    
            except RuntimeError:
                # Widgets have been deleted, ignore
                pass
                
        except Exception as e:
            print(f"Error updating output options state: {e}")
            # Ignore any other errors
            pass
            
    # ===== Deletion Canceler Options =====
    def _create_deletion_canceler_options(self):
        del_group = QGroupBox("Deletion Canceler Options")
        del_layout = QVBoxLayout()
        
        # Add the deletion canceler options to the container layout
        self.container_layout.addWidget(del_group)
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns where gaps may be due to sequencing artifacts "
            "rather than true biological deletions. These settings control how the cancellation is performed."
        )
        
        # Add the deletion canceler options to the container layout
        self.container_layout.addWidget(del_group)
        
        # Add stretch to push everything to the top
        self.container_layout.addStretch()
        
        # Set the container as the widget for the scroll area
        scroll.setWidget(container)
        
        # Add the scroll area to the page layout
        layout = QVBoxLayout(self)
        layout.addWidget(scroll)
        
        # Set the layout for the page
        self.setLayout(layout)
        
        # Mark widgets as initialized
        self.widgets_initialized = True
        
        # Connect signals last to avoid triggering updates during initialization
        if hasattr(self.wizard(), 'input_page'):
            if hasattr(self.wizard().input_page, 'species_phenotypes'):
                self.wizard().input_page.species_phenotypes.path_changed.connect(
                    self.update_output_options_state
                )
    
    def _create_deletion_canceler_options(self):
        """Create the deletion canceler options section."""
        del_group = QGroupBox("Deletion Canceler Options")
        del_layout = QVBoxLayout()
        
        # Add explanatory text
        explanation = QLabel(
            "The deletion canceler identifies and removes alignment columns where gaps may be due to sequencing artifacts "
            "rather than true biological deletions. These settings control how the cancellation is performed."
        )
        explanation.setWordWrap(True)
        del_layout.addWidget(explanation)
        
        # Add your deletion canceler options here
        
        del_group.setLayout(del_layout)
        self.container_layout.addWidget(del_group)
    
    def _update_grid_type_view(self):
        """Update the UI and config based on the selected grid type (linear or log)."""
        is_logspace = self.logspace_btn.isChecked()
        
        # Update the stacked widgets
        self.lambda1_stack.setCurrentIndex(1 if is_logspace else 0)
        self.lambda2_stack.setCurrentIndex(1 if is_logspace else 0)
        
        # Update the config based on the selected grid type
        if is_logspace:
            # When switching to logspace, use num_log_points value
            self.config.num_points = self.num_log_points.value()
        else:
            # When switching to linear, use lambda_step value
            self.config.num_points = self.lambda_step.value()
            
        # Update the config's grid type
        self.config.grid_type = 'log' if is_logspace else 'linear'
    
    def _update_deletion_canceler_state(self, state):
        """
        Update the state of the minimum aligned pairs control based on the 'Only cancel partner deletions' checkbox.
        
        Args:
            state: The state of the checkbox (2 for checked, 0 for unchecked)
        """
        is_checked = state == 2  # 2 is Qt.Checked
        self.min_pairs.setEnabled(is_checked)
        
        # Update the config
        self.config.cancel_only_partner = is_checked
        
        # If disabling, we don't need to include min_pairs in the config
        if not is_checked:
            # Set a default value that will be ignored by the CLI
            self.config.min_pairs = 0
    
    def _update_penalty_type(self, penalty_type):
        """
        Update the UI based on the selected group penalty type.
        
        Args:
            penalty_type: The selected penalty type ('median', 'std', 'sqrt', or 'linear')
        """
        # Update the config
        self.config.group_penalty_type = penalty_type
        
        # For std type, set initial and final values to be the same
        if penalty_type == 'std':
            current_value = self.initial_gp_value.value()
            self.final_gp_value.setValue(current_value)
            self.config.initial_gp_value = current_value
            self.config.final_gp_value = current_value
            
        # Show/hide the range settings based on the penalty type
        should_show = penalty_type not in ['median', 'std']
        if hasattr(self, 'penalty_range_group'):
            self.penalty_range_group.setVisible(should_show)
    
    def _update_phenotype_names_state(self):
        """Enable/disable the phenotype names section based on the selected output option."""
        is_gene_ranks_only = self.genes_only_btn.isChecked()
        
        # Update the config
        self.config.output_genes_only = is_gene_ranks_only
        
        # Enable/disable the phenotype name inputs
        self.pheno_label.setEnabled(not is_gene_ranks_only)
        self.pheno_name1.setEnabled(not is_gene_ranks_only)
        self.pheno_name2.setEnabled(not is_gene_ranks_only)
        
        # Update visual styling
        if is_gene_ranks_only:
            self.pheno_label.setStyleSheet("color: gray;")
            self.pheno_name1.setStyleSheet("color: gray; background-color: #f8f8f8;")
            self.pheno_name2.setStyleSheet("color: gray; background-color: #f8f8f8;")
        else:
            self.pheno_label.setStyleSheet("")
            self.pheno_name1.setStyleSheet("")
            self.pheno_name2.setStyleSheet("")