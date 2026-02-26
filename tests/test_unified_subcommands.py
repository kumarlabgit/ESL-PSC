import subprocess
from pathlib import Path

import pytest


def _rust_cli() -> Path:
    root = Path(__file__).resolve().parents[1]
    candidate = root / "esl_psc_rs" / "target" / "release" / "esl-psc"
    if not candidate.is_file():
        pytest.skip("release esl-psc binary not built at esl_psc_rs/target/release/esl-psc")
    return candidate


def test_unified_subcommand_help_smoke():
    rust_cli = _rust_cli()
    for subcommand in ("pairs", "plot", "site-counter"):
        result = subprocess.run(
            [str(rust_cli), subcommand, "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"{subcommand} --help failed: {result.stderr}"
        assert "usage:" in result.stdout.lower()


def test_unified_site_counter_matches_python_cli(tmp_path):
    rust_cli = _rust_cli()
    root = Path(__file__).resolve().parents[1]
    alignments_dir = root / "photosynthesis_alignments"
    groups_file = root / "photo_multi_species_groups.txt"
    tree_file = root / "photo_tree.nwk"

    py_out = tmp_path / "site_counter_python.csv"
    rs_out = tmp_path / "site_counter_unified.csv"

    py_result = subprocess.run(
        [
            "python3",
            "-m",
            "esl_psc_cli.fast_scan_cli",
            "--alignments_dir",
            str(alignments_dir),
            "--species_groups_file",
            str(groups_file),
            "--tree_file",
            str(tree_file),
            "--output_path",
            str(py_out),
        ],
        capture_output=True,
        text=True,
        cwd=root,
    )
    assert py_result.returncode == 0, py_result.stderr

    rs_result = subprocess.run(
        [
            str(rust_cli),
            "site-counter",
            "--alignments_dir",
            str(alignments_dir),
            "--species_groups_file",
            str(groups_file),
            "--tree_file",
            str(tree_file),
            "--output_path",
            str(rs_out),
        ],
        capture_output=True,
        text=True,
        cwd=root,
    )
    assert rs_result.returncode == 0, rs_result.stderr

    assert py_out.read_text(encoding="utf-8") == rs_out.read_text(encoding="utf-8")
