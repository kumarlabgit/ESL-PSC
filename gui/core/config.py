"""
Configuration model for ESL-PSC analysis.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
import os

@dataclass
class ESLConfig:
    # ─── Input files ────────────────────────────────────────────────────────────
    alignments_dir: str = ""
    species_groups_file: str = ""
    species_phenotypes_file: str = ""
    prediction_alignments_dir: str = ""
    limited_genes_file: str = ""
    response_dir: str = ""

    # ─── Hyper-parameters ───────────────────────────────────────────────────────
    initial_lambda1: float = 0.01
    final_lambda1: float = 0.99
    lambda_step: float = 0.1

    initial_lambda2: float = 0.01
    final_lambda2: float = 0.99
    lambda2_step: float = 0.1

    group_penalty_type: str = "median"
    initial_gp_value: float = 1.0
    final_gp_value: float = 6.0
    gp_step: float = 1.0

    grid_type: str = 'log'  # 'log' or 'linear'
    num_points: int = 20  # Number of points for log grid, or step size for linear grid

    # ─── Phenotype names ────────────────────────────────────────────────────────
    pheno_name1: str = "Convergent"
    pheno_name2: str = "Control"

    # ─── Deletion-canceller ─────────────────────────────────────────────────────
    nix_full_deletions: bool = False
    cancel_only_partner: bool = True
    min_pairs: int = 1

    # ─── Output options ─────────────────────────────────────────────────────────
    output_dir: str = os.path.join(os.getcwd(), "esl_psc_output")
    output_file_base_name: str = "esl_psc_results"
    keep_raw_output: bool = False
    show_selected_sites: bool = False
    no_genes_output: bool = False
    no_pred_output: bool = False

    # ─── Plot options ───────────────────────────────────────────────────────────
    no_sps_plot: bool = True
    make_sps_plot: bool = False
    make_sps_kde_plot: bool = False

    # ─── Multi-matrix / null models ─────────────────────────────────────────────
    top_rank_frac: float = 0.01
    make_null_models: bool = False
    make_pair_randomized_null_models: bool = False
    num_randomized_alignments: int = 10

    # ─── Helpers ────────────────────────────────────────────────────────────────
    def _flag(self, switch: bool, name: str) -> List[str]:
        return [name] if switch else []

    # ─── Public API ─────────────────────────────────────────────────────────────
    def _flag(self, b: bool, name: str) -> List[str]:
        return [name] if b else []

    def to_cli_args(self) -> List[str]:
        """Build a list of CLI flags that mirrors the current state."""
        a = []
        
        # Required parameters
        if not self.output_file_base_name:
            raise ValueError("missing output_file_base_name")
            
        # Input directories and files
        if self.alignments_dir:
            a += ["--alignments_dir", self.alignments_dir]
        if self.prediction_alignments_dir:
            a += ["--prediction_alignments_dir", self.prediction_alignments_dir]
        if self.species_groups_file:
            a += ["--species_groups_file", self.species_groups_file]
        if self.species_phenotypes_file:
            a += ["--species_pheno_path", self.species_phenotypes_file]
        if self.response_dir:
            a += ["--response_dir", self.response_dir]
            
        # Output configuration
        a += ["--output_file_base_name", self.output_file_base_name]
        
        # Sparsity parameters
        a += [
            "--initial_lambda1", str(self.initial_lambda1),
            "--final_lambda1", str(self.final_lambda1),
            "--initial_lambda2", str(self.initial_lambda2),
            "--final_lambda2", str(self.final_lambda2)
        ]
        
        if self.grid_type == 'log':
            a += ["--use_logspace", "--num_log_points", str(self.num_points)]
        else:  # linear
            a += ["--lambda_step", str(self.num_points)]  # Using num_points as step size for linear grid
            
        # Group penalty - only include values if not using 'median' or 'standard' type
        a += ["--group_penalty_type", self.group_penalty_type]
        if self.group_penalty_type not in ['median', 'standard']:
            a += [
                "--initial_gp_value", str(self.initial_gp_value),
                "--final_gp_value", str(self.final_gp_value),
                "--gp_step", str(self.gp_step)
            ]
        
        # Toggles
        a += self._flag(self.cancel_only_partner, "--cancel_only_partner")
        a += self._flag(self.nix_full_deletions, "--nix_full_deletions")
        
        # Output configuration (continued)
        a += ["--output_dir", self.output_dir]
        
        # Only include min_pairs if cancel_only_partner is True
        if hasattr(self, 'cancel_only_partner') and self.cancel_only_partner:
            a += ["--min_pairs", str(self.min_pairs)]
            
        # Add top_rank_frac if it exists
        if hasattr(self, 'top_rank_frac') and self.top_rank_frac is not None:
            a += ["--top_rank_frac", str(self.top_rank_frac)]
        
        # Output toggles
        a += self._flag(getattr(self, 'keep_raw_output', False), "--keep_raw_output")
        a += self._flag(getattr(self, 'show_selected_sites', False), "--show_selected_sites")
        a += self._flag(getattr(self, 'no_genes_output', False), "--no_genes_output")
        no_pred = getattr(self, 'no_pred_output', False)
        a += self._flag(no_pred, "--no_pred_output")
        
        # Only include plot flags if we're generating prediction output
        if not no_pred:
            a += self._flag(getattr(self, 'make_sps_plot', False), "--make_sps_plot")
            a += self._flag(getattr(self, 'make_sps_kde_plot', False), "--make_sps_kde_plot")
        
        # Null model options
        if hasattr(self, 'make_null_models') and self.make_null_models:
            a.append("--make_null_models")
            
        if hasattr(self, 'make_pair_randomized_null_models') and self.make_pair_randomized_null_models:
            a.append("--make_pair_randomized_null_models")
            if hasattr(self, 'num_randomized_alignments'):
                a += ["--num_randomized_alignments", str(self.num_randomized_alignments)]
        
        return a

    def get_command_string(self) -> str:
        """Return a full shell-ready command with proper quoting and formatting.
        
        Returns:
            str: The command string with each flag and its value on the same line,
                 and each flag-value pair on a new line.
        """
        parts = ["python -m esl_multimatrix"]
        args = self.to_cli_args()
        
        # Process arguments, grouping flags with their values
        i = 0
        while i < len(args):
            arg = args[i]
            # If this is a flag and there's a value after it, group them on one line
            if arg.startswith('--') and i + 1 < len(args) and not args[i+1].startswith('--'):
                parts.append(f"{arg} {args[i+1]}")
                i += 2
            # For standalone flags (no value), put them on their own line
            else:
                parts.append(arg)
                i += 1
        
        # Join with newlines and proper indentation
        return " \\\n  ".join(parts)