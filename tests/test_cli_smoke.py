# tests/test_cli_smoke.py

import subprocess
import sys
from pathlib import Path


def test_demo_run_smoke(tmp_path):
    """Runs the CLI on the demo data and checks that output files are created."""
    project_root = Path(__file__).parent.parent
    alignments_dir = project_root / "photosynthesis_alignments"
    species_groups_file = project_root / "photo_single_LC_matrix_species_groups.txt"
    species_pheno_file = project_root / "photo_species_phenotypes_full.txt"

    assert alignments_dir.exists(), "Demo alignments directory not found!"
    assert species_groups_file.exists(), "Demo species groups file not found!"

    output_basename = "smoke_test_output"
    output_dir = tmp_path

    command = [
        sys.executable,
        "-m", "esl_multimatrix",
        "--alignments_dir", str(alignments_dir),
        "--species_groups_file", str(species_groups_file),
        "--species_pheno_path", str(species_pheno_file),
        "--output_file_base_name", output_basename,
        "--output_dir", str(output_dir),
        "--use_logspace",
        "--num_log_points", "4",
        "--show_selected_sites",
    ]

    print(f"\nRunning command: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True)

    assert result.returncode == 0, (
        f"Script failed with exit code {result.returncode}\nStderr:\n{result.stderr}"
    )

    expected_predictions_file = output_dir / f"{output_basename}_species_predictions.csv"
    expected_ranks_file = output_dir / f"{output_basename}_gene_ranks.csv"
    expected_sites_file = output_dir / f"{output_basename}_selected_sites.csv"

    assert expected_predictions_file.exists(), "Predictions output file was not created."
    assert expected_ranks_file.exists(), "Gene ranks output file was not created."
    assert expected_sites_file.exists(), "Selected sites output file was not created."

    assert expected_ranks_file.stat().st_size > 0, "Gene ranks file is empty."
    assert expected_sites_file.stat().st_size > 0, "Selected sites file is empty."

    print("Smoke test passed: Script ran successfully and created output files.")
