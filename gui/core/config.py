"""
Configuration model for ESL-PSC analysis.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Any
import os

@dataclass
class ESLConfig:
    """Configuration class to store ESL-PSC analysis parameters."""
    
    # Input Files
    alignment_dir: str = ""
    species_groups_file: str = ""
    species_pheno_path: str = ""
    prediction_alignments_dir: str = ""
    limited_genes_list: str = ""
    canceled_alignments_dir: str = ""
    
    # Hyperparameters
    # Lambda (Sparsity) Parameters
    initial_lambda1: float = 0.01  # Position sparsity initial
    final_lambda1: float = 0.99    # Position sparsity final
    lambda1_step: float = 0.1      # Position sparsity step
    
    initial_lambda2: float = 0.01  # Group sparsity initial
    final_lambda2: float = 0.99    # Group sparsity final
    lambda2_step: float = 0.1      # Group sparsity step
    
    # Group Penalty Settings
    group_penalty_type: str = "median"  # Options: "default", "sqrt", "median", "linear"
    initial_gp_value: float = 1.0         # Initial group penalty constant
    final_gp_value: float = 6.0           # Final group penalty constant
    gp_step: float = 1.0                  # Group penalty step
    
    # Logspace Settings
    use_logspace: bool = True
    num_log_points: int = 20
    
    # Phenotype Names
    pheno_name1: str = "C4"  # Positive phenotype name
    pheno_name2: str = "C3"  # Negative phenotype name
    
    # Model Filtering
    min_genes: int = 0  # Minimum number of genes a model must have
    
    # Deletion Canceler Options
    nix_full_deletions: bool = False
    cancel_only_partner: bool = True
    min_pairs: int = 1
    
    # Output Options
    output_dir: str = os.path.join(os.getcwd(), 'esl_psc_output')
    output_file_base_name: str = "esl_psc_analysis"
    
    # Output Toggles
    keep_raw_output: bool = False
    show_selected_sites: bool = False
    no_genes_output: bool = False
    no_pred_output: bool = False
    
    # Plot Options
    make_sps_plot: bool = True
    make_sps_kde_plot: bool = False
    
    # Multi-matrix Options
    top_rank_frac: float = 0.01  # Fraction of top genes to highlight
    response_dir: str = ""      # Directory containing response matrices
    
    # Null Model Options
    make_null_models: bool = False
    make_pair_randomized_null_models: bool = False
    num_randomized_alignments: int = 10
    
    # Runtime Options
    num_threads: int = 1
    use_existing_preprocess: bool = False
    use_existing_alignments: bool = False
    use_uncanceled_alignments: bool = False
    delete_preprocess: bool = False
    
    # Additional Output Options
    save_model: bool = False
    save_predictions: bool = True
    generate_plots: bool = False
    plot_type: str = "violin"  # 'violin' or 'kde'
    
    def to_cli_args(self) -> List[str]:
        """Convert the configuration to a list of command-line arguments."""
        args = []
        
        # Required arguments
        if not self.alignment_dir:
            raise ValueError("Alignment directory is required")
        if not self.species_groups_file:
            raise ValueError("Species groups file is required")
            
        args.extend(["--alignment-dir", self.alignment_dir])
        args.extend(["--species-groups", self.species_groups_file])
        
        # Optional input files
        if self.species_phenotypes_file:
            args.extend(["--species-phenotypes", self.species_phenotypes_file])
        if self.prediction_alignments_dir:
            args.extend(["--prediction-alignments-dir", self.prediction_alignments_dir])
        if self.limited_genes_file:
            args.extend(["--limited-genes", self.limited_genes_file])
        
        # Parameters
        args.extend(["--lambda1", str(self.lambda1)])
        args.extend(["--lambda2", str(self.lambda2)])
        args.extend(["--group-penalty", self.group_penalty])
        args.extend(["--group-penalty-value", str(self.group_penalty_value)])
        args.extend(["--deletion-handling", self.deletion_handling])
        args.extend(["--min-species", str(self.min_species)])
        args.extend(["--min-aa", str(self.min_aa)])
        
        # Output options
        if self.output_dir:
            args.extend(["--output-dir", self.output_dir])
        if self.output_prefix != "esl_psc":
            args.extend(["--output-prefix", self.output_prefix])
        if self.save_model:
            args.append("--save-model")
        if self.save_predictions:
            args.append("--save-predictions")
        if self.generate_plots:
            args.append("--generate-plots")
            args.extend(["--plot-type", self.plot_type])
        
        return args
    
    def get_command_string(self) -> str:
        """Get the ESL-PSC command as a string."""
        return "python -m esl_multimatrix " + " ".join(f'"{arg}"' if ' ' in arg else arg 
                                                     for arg in self.to_cli_args())
