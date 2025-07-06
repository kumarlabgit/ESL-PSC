# tests/test_cli_smoke.py

import subprocess
import sys
from pathlib import Path


def test_demo_run_smoke(tmp_path):
    """Runs the CLI on the demo data and checks that output files are created."""
    project_root = Path(__file__).parent.parent
    alignments_dir = project_root / "photosynthesis_alignments"
    species_groups_src = project_root / "photo_single_LC_matrix_species_groups.txt"
    species_groups_file = tmp_path / "smoke_groups.txt"
    species_pheno_file = project_root / "photo_species_phenotypes.txt"

    assert alignments_dir.exists(), "Demo alignments directory not found!"
    assert species_groups_src.exists(), "Demo species groups file not found!"

    # Create a smaller species groups file using only the first two pairs
    with open(species_groups_src) as src, open(species_groups_file, "w") as dst:
        for _ in range(4):
            dst.write(src.readline())

    output_basename = "smoke_test_output"
    output_dir = tmp_path / "output"
    output_dir.mkdir()

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

    # Intermediate folders should be inside the output directory
    preprocess_dir = output_dir / "preprocessed_data_and_models"
    gap_dir = output_dir / f"{species_groups_file.stem}_gap-canceled_alignments"
    resp_dir = output_dir / f"{species_groups_file.stem}_response_matrices"

    for folder in (preprocess_dir, gap_dir, resp_dir):
        assert folder.exists() and any(folder.iterdir()), f"{folder} missing or empty"

    assert expected_ranks_file.stat().st_size > 0, "Gene ranks file is empty."
    assert expected_sites_file.stat().st_size > 0, "Selected sites file is empty."

    print("Smoke test passed: Script ran successfully and created output files.")
