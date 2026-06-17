from gui.core import fast_scan
from esl_psc_cli.esl_psc_functions import count_var_sites
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def write_fasta(path: str, records):
    with open(path, 'w') as f:
        for name, seq in records:
            f.write(f">{name}\n")
            f.write(seq + "\n")

def test_fast_scan_averages_fractional(tmp_path):
    # Create alignment directory
    align_dir = tmp_path / "alignments"
    align_dir.mkdir()

    # Create a single gene alignment with 6 species (A, B, C, D, E, OG)
    # Length 5. Designed so that:
    # - Combo1 (A,B vs C,D, outgroup OG) has 3 CCS sites
    # - Combo2 (A,E vs C,D, outgroup OG) has 2 CCS sites
    # Therefore average true CCS per gene should be (3 + 2) / 2 = 2.5
    records = [
        ("A",  "BBBAA"),
        ("B",  "BBBAA"),
        ("C",  "AAAAA"),
        ("D",  "AAAAA"),
        # E chosen so that (A,E) have B at positions 2 and 3 only (not at 1),
        # yielding 2 CCS for combo2 while combo1 has 3 CCS at positions 1..3.
        ("E",  "CBBAA"),
        ("OG", "AAAAA"),
    ]
    write_fasta(str(align_dir / "gene1.fas"), records)

    # Species groups with two combos:
    # 1) Convergent: A,B; Control: C,D
    # 2) Convergent: A,E; Control: C,D
    groups_path = tmp_path / "species_groups.txt"
    with open(groups_path, 'w') as g:
        g.write("A,B\n")
        g.write("C,D\n")
        g.write("A,E\n")
        g.write("C,D\n")

    # Run site counter
    results = fast_scan.fast_scan_alignments(str(align_dir), str(groups_path), "OG", progress_cb=None)
    assert results, "No results returned from fast_scan"
    # Find gene1
    rec = next((r for r in results if r.get("gene") == "gene1"), None)
    assert rec is not None, f"gene1 not found in results: {results}"

    # Expect fractional average true CCS: 2.5
    assert abs(rec["avg_true"] - 2.5) < 1e-6, f"Expected avg_true=2.5, got {rec['avg_true']}"
    # No control-convergence sites in this synthetic design
    assert abs(rec["avg_control"] - 0.0) < 1e-6, f"Expected avg_control=0.0, got {rec['avg_control']}"

    # Variable sites count should match CLI logic
    expected_var = count_var_sites([SeqRecord(Seq(seq), id=name) for name, seq in records])
    assert rec["variable_sites"] == expected_var
