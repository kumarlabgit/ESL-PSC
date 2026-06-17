import csv
import json
import subprocess
from pathlib import Path


def test_site_counter_cli_quotes_list_columns(tmp_path):
    root = Path(__file__).resolve().parents[1]
    alignments_dir = root / "demo_data" / "photosynthesis" / "alignments"
    groups_file = root / "demo_data" / "photosynthesis" / "photo_multi_species_groups.txt"
    tree_file = root / "demo_data" / "photosynthesis" / "photo_tree.nwk"
    out_csv = tmp_path / "site_counter.csv"

    result = subprocess.run(
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
            str(out_csv),
        ],
        capture_output=True,
        text=True,
        cwd=root,
    )
    assert result.returncode == 0, result.stderr

    with out_csv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    assert rows, "Site Counter CLI did not write any rows"
    first = rows[0]
    assert "per_combo_true" in first
    assert "per_combo_diff" in first
    assert "num_combos_top_frac" in first

    per_combo_true = json.loads(first["per_combo_true"])
    per_combo_diff = json.loads(first["per_combo_diff"])
    assert isinstance(per_combo_true, list)
    assert isinstance(per_combo_diff, list)
    assert len(per_combo_true) == len(per_combo_diff)

    num_top = float(first["num_combos_top_frac"])
    num_top_diff = float(first["num_combos_top_frac_by_diff"])
    num_top_ratio = float(first["num_combos_top_frac_by_ratio"])
    assert num_top >= 0
    assert num_top_diff >= 0
    assert num_top_ratio >= 0
