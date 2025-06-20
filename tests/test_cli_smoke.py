# tests/test_cli_smoke.py

import subprocess
import sys
from pathlib import Path

def test_demo_run_smoke(tmp_path):
    """
    Runs the full CLI script with the demo data and checks for success.
    This acts as a high-level integration smoke test.
    """
    # 1. Define paths relative to the project root
    project_root = Path(__file__).parent.parent
    alignments_dir = project_root / "photosynthesis_alignments"
    species_groups_file = project_root / "photo_single_LC_matrix_species_groups.txt"
    species_pheno_file = project_root / "photo_species_phenotypes.txt"
    
    # Ensure the required demo data exists before running the test
    assert alignments_dir.exists(), "Demo alignments directory not found!"
    assert species_groups_file.exists(), "Demo species groups file not found!"

    # 2. Construct the command arguments for the CLI script
    output_basename = "smoke_test_output"
    output_dir = tmp_path
    
    command = [
        sys.executable,  # Use the same python interpreter that's running the test
        "-m", "esl_multimatrix",
        "--alignments_dir", str(alignments_dir),
        "--species_groups_file", str(species_groups_file),
        "--species_pheno_path", str(species_pheno_file),
        "--output_file_base_name", output_basename,
        "--output_dir", str(output_dir),
        "--use_logspace",
        "--num_log_points", "4",  # Use a small number for a fast test
    ]

    # 3. Run the script as a subprocess
    print(f"\nRunning command: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True)

    # 4. Assert the results
    # Assert that the script completed successfully (exit code 0)
    assert result.returncode == 0, f"Script failed with exit code {result.returncode}\nStderr:\n{result.stderr}"

    # Assert that the expected output files were created
    expected_predictions_file = output_dir / f"{output_basename}_species_predictions.csv"
    expected_ranks_file = output_dir / f"{output_basename}_gene_ranks.csv"
    
    assert expected_predictions_file.exists(), "Predictions output file was not created."
    assert expected_ranks_file.exists(), "Gene ranks output file was not created."

    # Optional: A more detailed check to ensure the file is not empty
    assert expected_ranks_file.stat().st_size > 0, "Gene ranks file is empty."
    
    print("Smoke test passed: Script ran successfully and created output files.")