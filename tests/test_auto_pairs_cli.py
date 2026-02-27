import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("Bio")


def test_auto_pairs_cli_smoke(tmp_path):
    project_root = Path(__file__).parent.parent
    tree_file = project_root / "test_data" / "photosynthesis" / "photo_tree.nwk"
    pheno_file = project_root / "test_data" / "photosynthesis" / "photo_species_phenotypes.txt"
    out_file = tmp_path / "auto_pairs_groups.txt"

    assert tree_file.exists()
    assert pheno_file.exists()

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        str(out_file),
        "--method",
        "default",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"auto_pairs_cli failed: {result.stderr}"
    assert out_file.exists() and out_file.stat().st_size > 0

    lines = [ln.strip() for ln in out_file.read_text().splitlines() if ln.strip()]
    assert len(lines) % 2 == 0
    assert len(lines) >= 4


def test_auto_pairs_cli_num_random_sets(tmp_path):
    project_root = Path(__file__).parent.parent
    tree_file = project_root / "test_data" / "photosynthesis" / "photo_tree.nwk"
    pheno_file = project_root / "test_data" / "photosynthesis" / "photo_species_phenotypes.txt"
    out_base = tmp_path / "random_sets" / "pairs.txt"

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        str(out_base),
        "--method",
        "random",
        "--num_random_sets",
        "3",
        "--seed",
        "7",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"auto_pairs_cli random sets failed: {result.stderr}"

    expected = [
        tmp_path / "random_sets" / "pairs_001.txt",
        tmp_path / "random_sets" / "pairs_002.txt",
        tmp_path / "random_sets" / "pairs_003.txt",
    ]
    for path in expected:
        assert path.exists() and path.stat().st_size > 0


def test_auto_pairs_cli_num_random_sets_requires_random_method(tmp_path):
    project_root = Path(__file__).parent.parent
    tree_file = project_root / "test_data" / "photosynthesis" / "photo_tree.nwk"
    pheno_file = project_root / "test_data" / "photosynthesis" / "photo_species_phenotypes.txt"
    out_file = tmp_path / "pairs.txt"

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        str(out_file),
        "--method",
        "default",
        "--num_random_sets",
        "2",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    assert "--num_random_sets > 1 requires --method random" in result.stderr


def test_auto_pairs_cli_num_random_sets_output_directory(tmp_path):
    project_root = Path(__file__).parent.parent
    tree_file = project_root / "test_data" / "photosynthesis" / "photo_tree.nwk"
    pheno_file = project_root / "test_data" / "photosynthesis" / "photo_species_phenotypes.txt"
    out_dir = tmp_path / "random_sets_dir"

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        f"{out_dir}/",
        "--method",
        "random",
        "--num_random_sets",
        "2",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"auto_pairs_cli random dir sets failed: {result.stderr}"
    assert (out_dir / "auto_pairs_groups_001.txt").exists()
    assert (out_dir / "auto_pairs_groups_002.txt").exists()


def test_auto_pairs_cli_pct_contrast_smoke(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
    pheno_file = tmp_path / "pheno.csv"
    pheno_file.write_text(
        "species,value\nA,1\nB,2\nC,1\nD,2\n",
        encoding="utf-8",
    )
    out_file = tmp_path / "pct_groups.txt"

    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        str(out_file),
        "--method",
        "pct_contrast",
        "--min_pct_diff",
        "50",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"pct_contrast failed: {result.stderr}"
    lines = [ln.strip() for ln in out_file.read_text().splitlines() if ln.strip()]
    assert len(lines) >= 4 and len(lines) % 2 == 0


def test_auto_pairs_cli_pct_contrast_requires_positive_values(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
    pheno_file = tmp_path / "pheno.csv"
    pheno_file.write_text(
        "species,value\nA,0\nB,2\nC,1\nD,2\n",
        encoding="utf-8",
    )
    out_file = tmp_path / "pct_groups.txt"
    cmd = [
        sys.executable,
        "-m",
        "esl_psc_cli.auto_pairs_cli",
        "--tree_file",
        str(tree_file),
        "--species_pheno_path",
        str(pheno_file),
        "--output_path",
        str(out_file),
        "--method",
        "pct_contrast",
        "--min_pct_diff",
        "50",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    assert "strictly positive phenotype values" in result.stderr
