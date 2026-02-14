from pathlib import Path
from types import SimpleNamespace

from Bio import SeqIO

from esl_psc_cli import esl_multimatrix as emm
from esl_psc_cli import esl_psc_functions as ecf


def _write_multiline_fasta(path: Path) -> None:
    path.write_text(
        ">sp1\n"
        "ACGT\n"
        "AC\n"
        ">sp2\n"
        "TT\n"
        "GG\n",
        encoding="utf-8",
    )


def test_convert_alignment_dir_to_two_line_default_sibling(tmp_path):
    src = tmp_path / "alignments"
    src.mkdir()
    _write_multiline_fasta(src / "gene1.fas")

    out_dir, n_written = ecf.convert_alignment_dir_to_two_line(str(src))
    assert n_written == 1
    assert out_dir == str(tmp_path / "alignments_2line")
    out_file = Path(out_dir) / "gene1.fas"
    assert out_file.is_file()
    assert ecf.is_two_line_fasta(str(out_file))

    recs = list(SeqIO.parse(str(out_file), "fasta"))
    assert str(recs[0].seq) == "ACGTAC"
    assert str(recs[1].seq) == "TTGG"


def test_convert_alignment_dir_to_two_line_recursive_preserves_subdirs(tmp_path):
    src = tmp_path / "alns"
    nested = src / "combo_1"
    nested.mkdir(parents=True)
    _write_multiline_fasta(nested / "gene2.fasta")

    out_dir, n_written = ecf.convert_alignment_dir_to_two_line(
        str(src),
        recursive=True,
    )
    assert n_written == 1
    out_file = Path(out_dir) / "combo_1" / "gene2.fasta"
    assert out_file.is_file()
    assert ecf.is_two_line_fasta(str(out_file))


def test_find_non_two_line_fasta_files_detects_multiline(tmp_path):
    src = tmp_path / "alns"
    src.mkdir()
    _write_multiline_fasta(src / "gene3.fa")

    non_two = ecf.find_non_two_line_fasta_files(str(src))
    assert len(non_two) == 1
    assert non_two[0].endswith("gene3.fa")


def test_ensure_two_line_requires_auto_flag_when_multiline(tmp_path):
    src = tmp_path / "alignments"
    src.mkdir()
    _write_multiline_fasta(src / "gene4.fas")

    args = SimpleNamespace(auto_convert_to_2line=False)
    try:
        emm._ensure_two_line_with_optional_conversion(
            args,
            str(src),
            recursive=False,
            label="prediction alignments",
        )
        assert False, "Expected ValueError when auto-convert flag is missing"
    except ValueError as e:
        msg = str(e)
        assert "--auto_convert_to_2line" in msg
        assert "non-interactively" in msg


def test_ensure_two_line_errors_if_target_already_exists(tmp_path):
    src = tmp_path / "alignments"
    src.mkdir()
    _write_multiline_fasta(src / "gene5.fas")
    existing_target = tmp_path / "alignments_2line"
    existing_target.mkdir()

    args = SimpleNamespace(auto_convert_to_2line=True)
    try:
        emm._ensure_two_line_with_optional_conversion(
            args,
            str(src),
            recursive=False,
            label="prediction alignments",
        )
        assert False, "Expected ValueError when target _2line directory already exists"
    except ValueError as e:
        msg = str(e)
        assert "target directory already exists" in msg
        assert "Use that converted directory directly" in msg
