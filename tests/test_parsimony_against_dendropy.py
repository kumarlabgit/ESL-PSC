import pytest


pytest.importorskip("dendropy")

from Bio.Phylo.BaseTree import Clade, Tree
from Bio import Phylo

from gui.core.ancestral_reconstruction import (
    prune_tree,
    find_mrca,
    reconstruct_ancestral_sequence_with_sets,
)


def _write_newick(root: Clade, path: str) -> None:
    Phylo.write(Tree(root=root), path, "newick")


def _simple_outgroup_tree() -> Clade:
    # (((A,B)AB,(C,D)CD)ABCD,E)root
    A = Clade(name="A")
    B = Clade(name="B")
    C = Clade(name="C")
    D = Clade(name="D")
    E = Clade(name="E")

    AB = Clade(name="AB", clades=[A, B])
    CD = Clade(name="CD", clades=[C, D])
    ABCD = Clade(name="ABCD", clades=[AB, CD])
    return Clade(name="root", clades=[ABCD, E])


def test_mrca_state_sets_match_dendropy(tmp_path):
    # Small toy where some sites are ambiguous at the MRCA.
    tree_path = tmp_path / "t.nwk"
    _write_newick(_simple_outgroup_tree(), str(tree_path))

    seqs = {
        "A": "GATTA-A",
        "B": "GATTA-A",
        "C": "AATTA-A",
        "D": "AATTA-A",
        "E": "AAAAAAA",
    }
    analysis_species = {"A", "B", "C", "D"}

    # Our MRCA sets (downpass)
    full_tree = Phylo.read(str(tree_path), "newick").root
    pruned = prune_tree(full_tree, set(seqs.keys()))
    mrca = find_mrca(pruned, set(seqs.keys()) & analysis_species)
    _rep, ours_sets = reconstruct_ancestral_sequence_with_sets(pruned, seqs, mrca)

    # DendroPy MRCA sets (downpass reference)
    import dendropy
    from dendropy.model import parsimony

    labels = sorted(seqs.keys())
    taxa = dendropy.TaxonNamespace(labels)
    dtree = dendropy.Tree.get(
        path=str(tree_path),
        schema="newick",
        taxon_namespace=taxa,
        preserve_underscores=True,
        rooting="force-rooted",
    )
    dtree.retain_taxa_with_labels(labels)
    dtree.suppress_unifurcations()
    dtree.is_rooted = True

    alphabet = sorted({c for s in seqs.values() for c in s if c not in {"-", "?", "X", "x"}})
    code_of = {c: i for i, c in enumerate(alphabet)}
    char_of = {i: c for c, i in code_of.items()}
    seq_len = max(len(s) for s in seqs.values())

    taxon_state_sets_map = {}
    for taxon in taxa:
        s = seqs[taxon.label]
        ssl = []
        for pos in range(seq_len):
            ch = s[pos] if pos < len(s) else "?"
            if ch in {"-", "?", "X", "x"}:
                ssl.append(set())
            else:
                ssl.append({code_of[ch]})
        taxon_state_sets_map[taxon] = ssl

    parsimony.fitch_down_pass(
        dtree.postorder_node_iter(),
        state_sets_attr_name="ss",
        taxon_state_sets_map=taxon_state_sets_map,
    )
    dmrca = dtree.mrca(taxon_labels=sorted(set(labels) & set(analysis_species)))
    assert dmrca is not None
    ref_sets = [{char_of[i] for i in s} for s in getattr(dmrca, "ss")]

    assert len(ours_sets) == len(ref_sets)
    for a, b in zip(ours_sets, ref_sets):
        assert set(a) == set(b)
