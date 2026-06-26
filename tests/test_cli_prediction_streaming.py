import gzip
import json
import subprocess
from pathlib import Path

import pytest


def _rust_cli() -> Path:
    root = Path(__file__).resolve().parents[1]
    candidate = root / "esl_psc_rs" / "target" / "release" / "esl-psc"
    if not candidate.is_file():
        pytest.skip("release esl-psc binary not built at esl_psc_rs/target/release/esl-psc")
    return candidate


def _write_small_null_inputs(tmp_path: Path) -> tuple[Path, Path]:
    root = Path(__file__).resolve().parents[1]
    groups_src = root / "demo_data" / "photosynthesis" / "photo_single_LC_matrix_species_groups.txt"
    groups = tmp_path / "groups.txt"
    groups.write_text("".join(groups_src.read_text(encoding="utf-8").splitlines(True)[:8]))

    genes = tmp_path / "genes.txt"
    genes.write_text("ndhA.fas\nmatk.fas\n", encoding="utf-8")
    return groups, genes


def _run_small_response_flip_null(tmp_path: Path, output_dir: Path, base_name: str) -> subprocess.CompletedProcess:
    root = Path(__file__).resolve().parents[1]
    groups, genes = _write_small_null_inputs(tmp_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    return subprocess.run(
        [
            str(_rust_cli()),
            "--alignments_dir",
            str(root / "demo_data" / "photosynthesis" / "alignments"),
            "--species_groups_file",
            str(groups),
            "--species_pheno_path",
            str(root / "demo_data" / "photosynthesis" / "photo_species_phenotypes.txt"),
            "--limited_genes_list",
            str(genes),
            "--output_file_base_name",
            base_name,
            "--output_dir",
            str(output_dir),
            "--use_logspace",
            "--num_log_points",
            "2",
            "--make_null_models",
        ],
        capture_output=True,
        text=True,
        cwd=root,
    )


def test_checkpointed_predictions_are_sharded_and_not_serialized(tmp_path):
    output_dir = tmp_path / "out"
    result = _run_small_response_flip_null(tmp_path, output_dir, "streaming")
    assert result.returncode == 0, result.stderr

    predictions = output_dir / "streaming_species_predictions.csv"
    assert predictions.exists()
    assert sum(1 for _ in predictions.open(encoding="utf-8")) == 661

    checkpoint = output_dir / "checkpoint"
    shards = sorted(checkpoint.glob("predictions_combo_*.csv"))
    tmp_shards = sorted(checkpoint.glob("predictions_combo_*.csv.tmp"))
    assert len(shards) == 3
    assert tmp_shards == []

    with gzip.open(checkpoint / "state.json.gz", "rt", encoding="utf-8") as handle:
        state = json.load(handle)
    assert state["next_combo_index"] == 3
    assert state["prediction_rows"] == []


def test_prediction_shards_resume_without_duplicate_rows(tmp_path):
    output_dir = tmp_path / "out"
    result = _run_small_response_flip_null(tmp_path, output_dir, "resume")
    assert result.returncode == 0, result.stderr

    predictions = output_dir / "resume_species_predictions.csv"
    expected = predictions.read_text(encoding="utf-8")

    # Simulate a clean interruption after combo 0: combo 0's complete shard
    # remains, later shards/final CSV are absent, and meta.txt points to combo 0
    # as the last completed combo.
    for stale in (
        output_dir / "checkpoint" / "predictions_combo_1.csv",
        output_dir / "checkpoint" / "predictions_combo_2.csv",
        predictions,
    ):
        stale.unlink()
    (output_dir / "checkpoint" / "meta.txt").write_text("0\n", encoding="utf-8")

    resumed = _run_small_response_flip_null(tmp_path, output_dir, "resume")
    assert resumed.returncode == 0, resumed.stderr
    assert predictions.read_text(encoding="utf-8") == expected
    assert len(list((output_dir / "checkpoint").glob("predictions_combo_*.csv"))) == 3
