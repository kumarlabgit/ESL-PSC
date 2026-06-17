import subprocess
from pathlib import Path
from types import SimpleNamespace
from esl_psc_cli import deletion_canceler as dc
from esl_psc_cli import esl_psc_functions as ecf


def run_and_compare(tmp_path, extra_rs=None, extra_py=None):
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "tests" / "data" / "small_alignment"
    alignments_dir = data_dir / "alignments"
    response_file = data_dir / "response.txt"

    py_out = tmp_path / "python_out"
    rs_out = tmp_path / "rust_out"
    py_out.mkdir()
    rs_out.mkdir()

    species_list = ecf.get_species_to_check(str(response_file), check_order=True)
    args = SimpleNamespace(
        alignments_dir=str(alignments_dir),
        canceled_alignments_dir=str(py_out),
        nix_full_deletions=False,
        outgroup_species=None,
        cancel_only_partner=False,
        min_pairs=2,
        cancel_tri_allelic=False,
        limited_genes_list=None,
    )
    if extra_py:
        for k, v in extra_py.items():
            setattr(args, k, v)
    dc.generate_gap_canceled_alignments(args, [species_list])

    rs_cmd = [
        "cargo",
        "run",
        "--quiet",
        "--release",
        "--manifest-path",
        str(project_root / "legacy" / "deletion_canceler_rs" / "Cargo.toml"),
        "--",
        "--alignments-dir",
        str(alignments_dir),
        "--response-file",
        str(response_file),
        "--canceled-alignments-dir",
        str(rs_out),
    ]
    if extra_rs:
        rs_cmd.extend(extra_rs)
    subprocess.run(rs_cmd, check=True)

    py_file = py_out / "gene1.fas"
    rs_file = rs_out / "gene1.fas"
    assert py_file.read_text() == rs_file.read_text()


def test_rust_matches_python_default(tmp_path):
    run_and_compare(tmp_path)


def test_rust_matches_python_with_options(tmp_path):
    extra_rs = [
        "--cancel-only-partner",
        "--min-pairs",
        "1",
        "--cancel-tri-allelic",
        "--outgroup-species",
        "sp1",
    ]
    extra_py = {
        "cancel_only_partner": True,
        "min_pairs": 1,
        "cancel_tri_allelic": True,
        "outgroup_species": "sp1",
    }
    run_and_compare(tmp_path, extra_rs=extra_rs, extra_py=extra_py)
