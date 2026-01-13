import subprocess
import sys
from pathlib import Path

import pytest


def test_auto_pairs_cli_smoke(tmp_path):
    pytest.importorskip("Bio")

    project_root = Path(__file__).parent.parent
    tree_file = project_root / "photo_tree.nwk"
    pheno_file = project_root / "photo_species_phenotypes.txt"
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
