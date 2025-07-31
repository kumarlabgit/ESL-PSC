"""Lightweight checkpoint manager for ESL-PSC multicombo runs.

Files created inside ``<output_dir>/checkpoint``:
    * command.json  – JSON-serialized ``vars(args)`` of the first run
    * meta.txt      – last completed combo index (int)
    * state.pkl     – pickle of ``(gene_objects_dict, master_run_list)``
    * runs.jsonl    – append-only per-run audit records (one JSON per ESLRun)

A subsequent invocation will resume only if the current ``vars(args)`` exactly
matches the stored command. Otherwise the caller must supply the CLI flag
``--force_from_beginning`` which will delete the checkpoint and start fresh.
"""

from __future__ import annotations

import json
import os
import pickle
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, Optional, Tuple
from collections import defaultdict

__all__ = ["Checkpointer"]


class Checkpointer:
    """Utility class for saving and loading per-combo checkpoints."""

    # ------------------------------------------------------------------
    # Construction & file paths
    # ------------------------------------------------------------------
    def __init__(self, output_dir: str) -> None:
        self.cp_dir = os.path.join(output_dir, "checkpoint")
        os.makedirs(self.cp_dir, exist_ok=True)

        self.runs_file = os.path.join(self.cp_dir, "runs.jsonl")
        self.meta_file = os.path.join(self.cp_dir, "meta.txt")
        self.cmd_file = os.path.join(self.cp_dir, "command.json")
        self.state_file = os.path.join(self.cp_dir, "state.pkl")

    # ------------------------------------------------------------------
    # Query helpers
    # ------------------------------------------------------------------
    def has_checkpoint(self) -> bool:
        """Return True if any checkpoint metadata exists."""
        return os.path.exists(self.meta_file)

    def get_last_combo(self) -> Optional[int]:
        """Return the last completed combo index, or None if none."""
        if not os.path.exists(self.meta_file):
            return None
        try:
            with open(self.meta_file, "r") as fh:
                return int(fh.read().strip())
        except Exception:
            return None

    def load_command(self) -> Optional[Dict[str, Any]]:
        """Return the saved command dict, or None if absent."""
        if not os.path.exists(self.cmd_file):
            return None
        with open(self.cmd_file, "r") as fh:
            return json.load(fh)

    def is_same_command(self, current_cmd: Dict[str, Any]) -> bool:
        """Return True if *current_cmd* matches the stored command ignoring
        benign / auto-mutated args. Also stores a diff list on self.last_diff
        for debugging."""
        stored = self.load_command()
        if stored is None:
            self.last_diff = ["<no stored command>"]
            return False

        ignore = {
            # set/modified internally or non-critical toggles
            "esl_inputs_outputs_dir",
            "esl_main_dir",
            "canceled_alignments_dir",
            "use_existing_preprocess",
            "delete_preprocess",
            "use_uncanceled_alignments",
            "only_pos_gss",  # handled separately
            "lambda1_only",  # handled separately
            "make_sps_plot",  # plot flags can be toggled on resume
            "make_sps_kde_plot",
            "no_checkpoint",
            "force_from_beginning",
        }
        filtered_stored = {k: v for k, v in stored.items() if k not in ignore}
        filtered_current = {k: v for k, v in current_cmd.items() if k not in ignore}

        def _values_equal(a, b, key):
            # Treat None and False as equivalent for specific flags
            if key in {"only_pos_gss", "lambda1_only"}:
                if (a is None and b is False) or (a is False and b is None):
                    return True
            # For alignments_dir, allow stored None/'' to equal any derived path
            if key == "alignments_dir":
                if a in {None, ""} or b in {None, ""}:
                    return True
            if key == "use_existing_alignments":
                # Allow auto-upgrade from False/None to True on resume
                if a in {None, False} and b is True:
                    return True
                if b in {None, False} and a is True:
                    return True
            return a == b

        # Keep a diff for diagnostics
        self.last_diff = [
            f"{k}: stored={filtered_stored.get(k)} current={filtered_current.get(k)}"
            for k in sorted(set(filtered_stored) | set(filtered_current))
            if not _values_equal(filtered_stored.get(k), filtered_current.get(k), k)
        ]
        return not self.last_diff

    def load_state(self) -> Tuple[Any, List[Any]]:
        """Load and return (gene_objects_dict, master_run_list)."""
        if not os.path.exists(self.state_file):
            return None, []
        with open(self.state_file, "rb") as fh:
            state = pickle.load(fh)

        gene_objects = state.get("gene_objects_dict")
        master_runs = state.get("master_run_list")

        # Restore defaultdict behavior stripped during checkpointing
        if gene_objects:
            for gene in gene_objects.values():
                if hasattr(gene, "selected_sites") and isinstance(gene.selected_sites, dict):
                    gene.selected_sites = defaultdict(lambda: 0, gene.selected_sites)

        if master_runs:
            for run in master_runs:
                if hasattr(run, "species_scores") and isinstance(run.species_scores, dict):
                    run.species_scores = defaultdict(lambda: 0.0, run.species_scores)

        return gene_objects, master_runs

    # ------------------------------------------------------------------
    # Writing helpers
    # ------------------------------------------------------------------
    def _append_run_audit(self, combo_idx: int, runs: List[Any]) -> None:
        if not runs:
            return
        with open(self.runs_file, "a") as fh:
            for run in runs:
                fh.write(
                    json.dumps(
                        {
                            "combo": combo_idx,
                            "lambda1": run.lambda1,
                            "lambda2": run.lambda2,
                            "penalty_term": run.run_family.penalty_term,
                            "input_rmse": run.input_rmse,
                        }
                    )
                    + "\n"
                )

    def save_checkpoint(
        self,
        combo_idx: int,
        gene_objects_dict: Any,
        master_run_list: List[Any],
        args_dict: Dict[str, Any],
        runs_this_combo: List[Any],
    ) -> None:
        """Persist full checkpoint after finishing *combo_idx*.

        Ensures the checkpoint directory exists (it may have been deleted by
        --force_from_beginning earlier)."""
        # Guarantee checkpoint directory exists (may have been removed)
        os.makedirs(self.cp_dir, exist_ok=True)
        # 1) append audit runs
        self._append_run_audit(combo_idx, runs_this_combo)

        # 2) atomically write pickle state
        # Convert any un-picklable defaultdict(lambda: 0) to a plain dict
        # during pickling but restore afterwards so in-memory objects retain
        # their convenient default behaviour.
        modified_genes = []
        for gene in gene_objects_dict.values():
            if hasattr(gene, "selected_sites") and isinstance(gene.selected_sites, defaultdict):
                gene.selected_sites = dict(gene.selected_sites)  # strip the lambda
                modified_genes.append(gene)
        # Also sanitise any defaultdicts in ESLRun objects
        modified_runs = []
        for run in master_run_list:
            if hasattr(run, "species_scores") and isinstance(run.species_scores, defaultdict):
                run.species_scores = dict(run.species_scores)
                modified_runs.append(run)
        tmp = NamedTemporaryFile("wb", delete=False, dir=self.cp_dir)
        pickle.dump(
            {
                "gene_objects_dict": gene_objects_dict,
                "master_run_list": master_run_list,
            },
            tmp,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
        tmp.flush()
        os.fsync(tmp.fileno())
        tmp.close()
        os.replace(tmp.name, self.state_file)

        # Restore defaultdict wrappers now that pickling is done
        for gene in modified_genes:
            gene.selected_sites = defaultdict(lambda: 0, gene.selected_sites)
        for run in modified_runs:
            run.species_scores = defaultdict(lambda: 0.0, run.species_scores)

        # 3) update last combo index
        with open(self.meta_file, "w") as fh:
            fh.write(str(combo_idx))

        # 4) save command once
        if not os.path.exists(self.cmd_file):
            with open(self.cmd_file, "w") as fh:
                json.dump(args_dict, fh, indent=2)

    # ------------------------------------------------------------------
    # Maintenance helpers
    # ------------------------------------------------------------------
    def clear(self) -> None:
        """Delete the entire checkpoint directory."""
        if not os.path.isdir(self.cp_dir):
            return
        for name in os.listdir(self.cp_dir):
            try:
                os.remove(os.path.join(self.cp_dir, name))
            except Exception:
                pass
        try:
            os.rmdir(self.cp_dir)
        except Exception:
            pass

    # Backwards-compat alias (old code called checkpoint_runs)
    def checkpoint_runs(
        self,
        combo_idx: int,
        runs: List[Any],
        gene_objects_dict: Any = None,
        master_run_list: List[Any] = None,
        args_dict: Dict[str, Any] = None,
    ) -> None:  # noqa: D401
        """Legacy wrapper – delegates to save_checkpoint with minimal params."""
        # Only writes audit; full save_checkpoint should be used by new code.
        self._append_run_audit(combo_idx, runs)

