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

    use_logspace: bool = True
    num_log_points: int = 20

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
        if self.response_dir:
            a += ["--response_dir", self.response_dir]
            
        # Output configuration
        a += ["--output_file_base_name", self.output_file_base_name]
        
        # Sparsity parameters
        if self.use_logspace:
            a += [
                "--use_logspace",
                "--initial_lambda1", str(self.initial_lambda1),
                "--final_lambda1", str(self.final_lambda1),
                "--initial_lambda2", str(self.initial_lambda2),
                "--final_lambda2", str(self.final_lambda2),
                "--num_log_points", str(self.num_log_points)
            ]
        else:
            a += [
                "--initial_lambda1", str(self.initial_lambda1),
                "--final_lambda1", str(self.final_lambda1),
                "--initial_lambda2", str(self.initial_lambda2),
                "--final_lambda2", str(self.final_lambda2),
                "--lambda_step", str(self.lambda_step)
            ]
            
        # Group penalty
        a += [
            "--group_penalty_type", self.group_penalty_type,
            "--initial_gp_value", str(self.initial_gp_value),
            "--final_gp_value", str(self.final_gp_value),
            "--gp_step", str(self.gp_step)
        ]
        
        # Toggles
        a += self._flag(self.cancel_only_partner, "--cancel_only_partner")
        a += self._flag(self.nix_full_deletions, "--nix_full_deletions")
        
        # Output configuration (continued)
        a += [
            "--min_pairs", str(self.min_pairs),
            "--output_dir", self.output_dir
        ]
        
        # Output toggles
        a += self._flag(self.keep_raw_output, "--keep_raw_output")
        a += self._flag(self.no_genes_output, "--no_genes_output")
        a += self._flag(self.no_pred_output, "--no_pred_output")
        a += self._flag(self.make_sps_plot, "--make_sps_plot")
        a += self._flag(self.make_sps_kde_plot, "--make_sps_kde_plot")
        
        return a

    def get_command_string(self) -> str:
        """Return a full shell-ready command with proper quoting and formatting.
        
        Returns:
            str: The command string with flag-value pairs on the same line.
        """
        parts = ["python", "-m", "esl_multimatrix"]
        args = self.to_cli_args()
        
        # Group flags with their values
        i = 0
        while i < len(args):
            part = args[i]
            # If this is a flag and there's a value after it, group them
            if part.startswith('--') and i + 1 < len(args) and not args[i+1].startswith('--'):
                parts.append(f"{part} {args[i+1]}")
                i += 2
            else:
                parts.append(part)
                i += 1
        
        # Format with line continuations
        cmd_lines = []
        current_line = []
        line_length = 0
        max_line_length = 80  # Target line length before wrapping
        
        for part in parts:
            # For the first part or if adding this part would exceed max line length
            if not current_line or line_length + len(part) + 1 > max_line_length:
                if current_line:  # If there's a current line, add it to the result
                    cmd_lines.append(" ".join(current_line))
                current_line = [part]
                line_length = len(part)
            else:
                current_line.append(part)
                line_length += len(part) + 1  # +1 for the space
        
        # Add the last line
        if current_line:
            cmd_lines.append(" ".join(current_line))
        
        # Join with line continuations and proper indentation
        return " \\\n  ".join(cmd_lines)