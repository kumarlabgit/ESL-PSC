import pytest

pytest.importorskip("Bio")

from esl_psc_cli.auto_pairs import load_tree, sweep_pct_contrast_pair_counts


def test_pct_contrast_sweep_counts_monotonic(tmp_path):
    tree_path = tmp_path / "tree.nwk"
    tree_path.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
    tree = load_tree(str(tree_path))
    phenos = {"A": 1.0, "B": 2.0, "C": 1.0, "D": 2.0}
    sweep = sweep_pct_contrast_pair_counts(tree, phenos, [0.0, 25.0, 50.0, 150.0])
    counts = [c for _t, c in sweep]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] >= 2
    assert counts[-1] == 0
