"""
Configuration model for ESL-PSC analysis.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Any

@dataclass
class ESLConfig:
    """Configuration for an ESL-PSC analysis."""
    # Input files
    alignment_dir: str = ""
    species_groups_file: str = ""
    species_phenotypes_file: str = ""
    prediction_alignments_dir: str = ""
    limited_genes_file: str = ""
    
    # Parameters
    lambda1: float = 0.1
    lambda2: float = 0.1
    group_penalty: str = "l1"  # 'l1', 'l2', or 'l1l2'
    group_penalty_value: float = 1.0
    deletion_handling: str = "drop"  # 'drop', 'impute', or 'impute_consensus'
    min_species: int = 1
    min_aa: int = 1
    
    # Output options
    output_dir: str = ""
    output_prefix: str = "esl_psc"
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
