import subprocess
import sys
from pathlib import Path

def test_continuous_pheno_disables_plots(tmp_path):
    project_root = Path(__file__).parent.parent
    alignments_dir = project_root / "photosynthesis_alignments"
    species_groups_src = project_root / "photo_single_LC_matrix_species_groups.txt"
    species_groups_file = tmp_path / "groups.txt"
    species_pheno_src = project_root / "photo_species_phenotypes.txt"
    continuous_pheno_file = tmp_path / "continuous_pheno.txt"
    limited_genes_file = tmp_path / "limited_genes.txt"

    assert alignments_dir.exists()
    assert species_groups_src.exists()

    # use first four species entries to keep run small
    with open(species_groups_src) as src, open(species_groups_file, "w") as dst:
        for _ in range(4):
            dst.write(src.readline())

    # top two genes for fast run
    top_genes = ["ndhA.fas", "matk.fas"]
    limited_genes_file.write_text("\n".join(top_genes))

    # build a continuous phenotype file using values 1.0..n*1.0
    cont_lines = []
    expected_values = {}
    with open(species_pheno_src) as src:
        for idx, line in enumerate(src, 1):
            line = line.strip()
            if not line:
                continue
            species = line.split(",")[0]
            value = float(idx)
            cont_lines.append(f"{species},{value}\n")
            expected_values[species] = value
    continuous_pheno_file.write_text("".join(cont_lines))

    output_basename = "continuous_test_output"
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cmd = [
        sys.executable,
        "-m", "esl_multimatrix",
        "--alignments_dir", str(alignments_dir),
        "--species_groups_file", str(species_groups_file),
        "--species_pheno_path", str(continuous_pheno_file),
        "--limited_genes_list", str(limited_genes_file),
        "--output_file_base_name", output_basename,
        "--output_dir", str(output_dir),
        "--use_continuous_phenotypes",
        "--use_logspace",
        "--num_log_points", "2",
        "--make_sps_plot",
        "--show_selected_sites",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "continuous response variables" in result.stdout

    # SPS plots should be skipped even if requested
    plot_file = output_dir / f"{output_basename}_pred_sps_plot.svg"
    assert not plot_file.exists()

    # Core output tables should still be produced
    predictions = output_dir / f"{output_basename}_species_predictions.csv"
    ranks = output_dir / f"{output_basename}_gene_ranks.csv"
    sites = output_dir / f"{output_basename}_selected_sites.csv"
    for f in (predictions, ranks, sites):
        assert f.exists() and f.stat().st_size > 0

    # Response matrix should contain the fabricated continuous values
    resp_dir = output_dir / f"{species_groups_file.stem}_response_matrices"
    resp_file = resp_dir / "combo_0.txt"
    assert resp_file.exists()
    contents = resp_file.read_text().strip().splitlines()
    for line in contents:
        species, value = line.split("\t")
        assert float(value) == expected_values[species]

    # Now ensure the continuous plot option works
    output_dir2 = tmp_path / "out2"
    output_dir2.mkdir()
    base2 = "continuous_plot_output"
    cmd2 = [
        sys.executable,
        "-m", "esl_multimatrix",
        "--alignments_dir", str(alignments_dir),
        "--species_groups_file", str(species_groups_file),
        "--species_pheno_path", str(continuous_pheno_file),
        "--limited_genes_list", str(limited_genes_file),
        "--output_file_base_name", base2,
        "--output_dir", str(output_dir2),
        "--use_continuous_phenotypes",
        "--use_logspace",
        "--num_log_points", "2",
        "--make_continuous_plot",
    ]
    result2 = subprocess.run(cmd2, capture_output=True, text=True)
    assert result2.returncode == 0, result2.stderr
    cont_plot = output_dir2 / f"{base2}_continuous_plot.svg"
    assert cont_plot.exists()


def test_response_dir_with_continuous_values(tmp_path):
    project_root = Path(__file__).parent.parent
    alignments_dir = project_root / "photosynthesis_alignments"
    species_pheno_src = project_root / "photo_species_phenotypes.txt"
    response_dir = tmp_path / "resp"
    response_dir.mkdir()

    species = []
    with open(species_pheno_src) as src:
        for line in src:
            if len(species) >= 4:
                break
            species.append(line.split(",")[0])
    with open(response_dir / "resp.txt", "w") as dst:
        for idx, sp in enumerate(species, 1):
            dst.write(f"{sp}\t{float(idx)}\n")

    limited_genes_file = tmp_path / "limited_genes.txt"
    limited_genes_file.write_text("ndhA.fas\nmatk.fas")

    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cmd = [
        sys.executable,
        "-m", "esl_multimatrix",
        "--alignments_dir", str(alignments_dir),
        "--response_dir", str(response_dir),
        "--limited_genes_list", str(limited_genes_file),
        "--output_file_base_name", "out",
        "--output_dir", str(output_dir),
        "--use_logspace",
        "--num_log_points", "2",
        "--no_pred_output",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "continuous response variables" in result.stdout
