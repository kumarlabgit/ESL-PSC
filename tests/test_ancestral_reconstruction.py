"""Tests for parsimony-based ancestral reconstruction for Fast Scan."""
import pytest
import tempfile
import os
from io import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from gui.core.ancestral_reconstruction import (
    parse_tree,
    get_terminal_names,
    prune_tree,
    find_mrca,
    is_mrca_at_root,
    fitch_parsimony_downpass,
    fitch_parsimony_uppass,
    reconstruct_ancestral_sequence,
    reconstruct_ancestral_sequence_with_sets,
    get_ancestral_outgroup_for_alignment,
    validate_tree_for_fast_scan,
    AncestralReconstructionError,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures and helper functions
# ─────────────────────────────────────────────────────────────────────────────

def create_simple_tree():
    """Create a simple tree: ((A,B)AB,(C,D)CD)root."""
    # Build clades
    A = Clade(name="A")
    B = Clade(name="B")
    C = Clade(name="C")
    D = Clade(name="D")
    
    AB = Clade(name="AB", clades=[A, B])
    CD = Clade(name="CD", clades=[C, D])
    
    root = Clade(name="root", clades=[AB, CD])
    
    return root


def create_tree_with_outgroup():
    """Create tree with outgroup: (((A,B)AB,(C,D)CD)ABCD,E)root."""
    A = Clade(name="A")
    B = Clade(name="B")
    C = Clade(name="C")
    D = Clade(name="D")
    E = Clade(name="E")
    
    AB = Clade(name="AB", clades=[A, B])
    CD = Clade(name="CD", clades=[C, D])
    ABCD = Clade(name="ABCD", clades=[AB, CD])
    
    root = Clade(name="root", clades=[ABCD, E])
    
    return root


def write_newick_to_file(tree: Clade, path: str):
    """Write tree to Newick file."""
    temp_tree = Tree(root=tree)
    Phylo.write(temp_tree, path, 'newick')


@pytest.fixture
def simple_tree_file(tmp_path):
    """Fixture providing a simple tree file."""
    tree = create_simple_tree()
    path = tmp_path / "simple.nwk"
    write_newick_to_file(tree, str(path))
    return str(path)


@pytest.fixture
def outgroup_tree_file(tmp_path):
    """Fixture providing a tree with outgroup."""
    tree = create_tree_with_outgroup()
    path = tmp_path / "outgroup.nwk"
    write_newick_to_file(tree, str(path))
    return str(path)


@pytest.fixture
def species_groups_file(tmp_path):
    """Create a simple species groups file."""
    path = tmp_path / "groups.txt"
    with open(path, 'w') as f:
        f.write("A\n")
        f.write("C\n")
        f.write("B\n")
        f.write("D\n")
    return str(path)


# ─────────────────────────────────────────────────────────────────────────────
# Test tree parsing
# ─────────────────────────────────────────────────────────────────────────────

def test_parse_tree_newick(simple_tree_file):
    """Test parsing a Newick tree file."""
    tree = parse_tree(simple_tree_file)
    assert tree is not None
    terminals = get_terminal_names(tree)
    assert "A" in terminals
    assert "B" in terminals
    assert "C" in terminals
    assert "D" in terminals


def test_parse_tree_nonexistent():
    """Test parsing a nonexistent file raises error."""
    with pytest.raises(AncestralReconstructionError, match="not found"):
        parse_tree("/nonexistent/path/tree.nwk")


def test_parse_tree_invalid():
    """Test parsing an invalid tree raises error."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write("this is not a tree")
        path = f.name
    
    try:
        with pytest.raises(AncestralReconstructionError, match="Failed to parse"):
            parse_tree(path)
    finally:
        os.unlink(path)


# ─────────────────────────────────────────────────────────────────────────────
# Test tree pruning
# ─────────────────────────────────────────────────────────────────────────────

def test_prune_tree_keep_all():
    """Test pruning tree to keep all species."""
    tree = create_simple_tree()
    keep = {"A", "B", "C", "D"}
    pruned = prune_tree(tree, keep)
    
    terminals = get_terminal_names(pruned)
    assert terminals == keep


def test_prune_tree_keep_subset():
    """Test pruning tree to keep subset of species."""
    tree = create_simple_tree()
    keep = {"A", "B"}
    pruned = prune_tree(tree, keep)
    
    terminals = get_terminal_names(pruned)
    assert terminals == keep
    assert "C" not in terminals
    assert "D" not in terminals


def test_prune_tree_collapse_single_child():
    """Test that single-child nodes are collapsed after pruning."""
    tree = create_tree_with_outgroup()
    keep = {"A", "E"}
    pruned = prune_tree(tree, keep)
    
    terminals = get_terminal_names(pruned)
    assert terminals == keep
    
    # Check that internal structure is simplified
    # The pruned tree should have minimal depth
    assert len(list(pruned.find_clades())) <= 3  # root + 2 terminals


def test_prune_tree_no_species_match():
    """Test pruning with no matching species raises error."""
    tree = create_simple_tree()
    keep = {"X", "Y", "Z"}
    
    with pytest.raises(AncestralReconstructionError, match="None of the requested species"):
        prune_tree(tree, keep)


# ─────────────────────────────────────────────────────────────────────────────
# Test MRCA identification
# ─────────────────────────────────────────────────────────────────────────────

def test_find_mrca_two_species():
    """Test finding MRCA of two species."""
    tree = create_simple_tree()
    mrca = find_mrca(tree, {"A", "B"})
    
    # MRCA should be the AB node
    descendants = get_terminal_names(mrca)
    assert "A" in descendants
    assert "B" in descendants
    assert "C" not in descendants
    assert "D" not in descendants


def test_find_mrca_all_species():
    """Test finding MRCA of all species (should be root)."""
    tree = create_simple_tree()
    mrca = find_mrca(tree, {"A", "B", "C", "D"})
    
    # MRCA should be the root
    assert mrca == tree


def test_find_mrca_single_species():
    """Test finding MRCA of single species returns that terminal."""
    tree = create_simple_tree()
    mrca = find_mrca(tree, {"A"})
    
    # MRCA should be A itself
    assert mrca.name == "A"
    assert mrca.is_terminal()


def test_find_mrca_no_species():
    """Test finding MRCA with empty set raises error."""
    tree = create_simple_tree()
    
    with pytest.raises(AncestralReconstructionError, match="No species provided"):
        find_mrca(tree, set())


def test_find_mrca_nonexistent_species():
    """Test finding MRCA with nonexistent species raises error."""
    tree = create_simple_tree()
    
    with pytest.raises(AncestralReconstructionError, match="No matching species"):
        find_mrca(tree, {"X", "Y"})


# ─────────────────────────────────────────────────────────────────────────────
# Test MRCA at root detection
# ─────────────────────────────────────────────────────────────────────────────

def test_is_mrca_at_root_true():
    """Test detecting when MRCA is at root."""
    tree = create_simple_tree()
    mrca = find_mrca(tree, {"A", "B", "C", "D"})
    
    assert is_mrca_at_root(tree, mrca) is True


def test_is_mrca_at_root_false():
    """Test detecting when MRCA is not at root."""
    tree = create_simple_tree()
    mrca = find_mrca(tree, {"A", "B"})
    
    assert is_mrca_at_root(tree, mrca) is False


def test_is_mrca_at_root_with_single_child_wrapper():
    """Test MRCA detection unwraps single-child root nodes."""
    # Create tree with single-child wrapper: (((A,B)AB,C)ABC)root
    A = Clade(name="A")
    B = Clade(name="B")
    C = Clade(name="C")
    AB = Clade(name="AB", clades=[A, B])
    ABC = Clade(name="ABC", clades=[AB, C])
    root = Clade(name="root", clades=[ABC])  # Single child
    
    mrca = find_mrca(root, {"A", "B", "C"})
    
    # Even though tree has a wrapper, MRCA is at effective root
    assert is_mrca_at_root(root, mrca) is True


def test_is_mrca_at_root_with_outgroup():
    """Test MRCA not at root when outgroup present."""
    tree = create_tree_with_outgroup()
    mrca = find_mrca(tree, {"A", "B", "C", "D"})
    
    # MRCA is ABCD clade, not root (E is outgroup)
    assert is_mrca_at_root(tree, mrca) is False


# ─────────────────────────────────────────────────────────────────────────────
# Test Fitch parsimony
# ─────────────────────────────────────────────────────────────────────────────

def test_fitch_parsimony_simple():
    """Test Fitch parsimony on a simple case."""
    tree = create_simple_tree()
    
    # Sequences where A and B share 'G', C and D share 'A'
    sequences = {
        "A": "G",
        "B": "G",
        "C": "A",
        "D": "A",
    }
    
    # Reconstruct at root
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    # Root should be either G or A (both are parsimonious)
    assert len(ancestral_seq) == 1
    assert ancestral_seq[0] in {"G", "A"}


def test_fitch_parsimony_unambiguous():
    """Test Fitch parsimony with unambiguous reconstruction."""
    tree = create_simple_tree()
    
    # All same residue
    sequences = {
        "A": "G",
        "B": "G",
        "C": "G",
        "D": "G",
    }
    
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    # Root should definitely be G
    assert ancestral_seq == "G"


def test_fitch_parsimony_with_gaps():
    """Test Fitch parsimony correctly handles gaps."""
    tree = create_simple_tree()
    
    # Sequences with gaps (should be ignored in parsimony)
    sequences = {
        "A": "G-A",
        "B": "G-T",
        "C": "-GA",
        "D": "TGT",
    }
    
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    # Should reconstruct 3 positions
    assert len(ancestral_seq) == 3
    # Position 0: G is common to A,B; T in D
    # Position 1: Gap in A,B; G in C,D
    # Position 2: A in A,C; T in B,D


def test_fitch_parsimony_missing_data():
    """Test Fitch parsimony with missing/ambiguous data."""
    tree = create_simple_tree()
    
    sequences = {
        "A": "G?X",
        "B": "GTX",
        "C": "ATC",
        "D": "ATC",
    }
    
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    assert len(ancestral_seq) == 3


def test_fitch_parsimony_longer_sequence():
    """Test Fitch parsimony on a longer sequence."""
    tree = create_simple_tree()
    
    sequences = {
        "A": "GATTACA",
        "B": "GATTACA",
        "C": "GATTAGA",
        "D": "GATTAGA",
    }
    
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    assert len(ancestral_seq) == 7
    # First 5 positions are identical across all species
    assert ancestral_seq[:5] == "GATTA"
    # Position 5: C in A,B; G in C,D (ambiguous)
    assert ancestral_seq[5] in {"C", "G"}
    # Position 6: A everywhere
    assert ancestral_seq[6] == "A"


# ─────────────────────────────────────────────────────────────────────────────
# New tests: MRCA state-set behavior and X/x handling
# ─────────────────────────────────────────────────────────────────────────────

def test_mrca_sets_ambiguous_no_count_under_tolerant():
    """When MRCA set is ambiguous and includes both control and convergent residues,
    tolerant set-based logic must NOT count the site as convergent.

    Tree: ((A,B),(C,D))root
    Sequences at position 0: A=B=G, C=D=A -> MRCA set {A,G}
    Controls=C,D=A (in MRCA set), Convergents=A,B=G (also in MRCA set) -> no count.
    """
    tree = create_simple_tree()  # ((A,B),(C,D))root
    sequences = {"A":"G", "B":"G", "C":"A", "D":"A"}
    # Reconstruct with sets
    # Wrap to full sequences of length 1
    rep, sets = reconstruct_ancestral_sequence_with_sets(tree, sequences, tree)
    assert len(sets) == 1
    mrca_set = sets[0]
    assert mrca_set == {"A", "G"}
    # Emulate tolerant gating
    ctrl = ["A","A"]
    conv = ["G","G"]
    min_agree = 1.0
    agree = sum(1 for r in ctrl if r in mrca_set) / len(ctrl)
    assert agree >= min_agree
    # Convergent derived must be NOT in MRCA set to count; here it is in the set
    assert all(r in mrca_set for r in conv)
    counted = any((r not in mrca_set) for r in conv) and conv.count(conv[0]) >= 2
    assert not counted


def test_mrca_sets_unambiguous_counts_when_derived_not_in_set():
    """When MRCA is unambiguous (e.g., {A}) and convergents share G, count CCS.

    Tree: ((A,B),(C,D))root
    Sequences: A=B=G; C=D=A -> MRCA set {A}; controls match A; convergents share G not in set -> count.
    """
    tree = create_simple_tree()
    sequences = {"A":"G", "B":"G", "C":"A", "D":"A"}
    rep, sets = reconstruct_ancestral_sequence_with_sets(tree, sequences, tree)
    mrca_set = sets[0]
    # Note: In this topology, the MRCA set is actually {A,G};
    # simulate an unambiguous scenario by changing D to A and C to A and add an outgroup E=A to drive parsimony to A.
    from Bio.Phylo.BaseTree import Clade
    E = Clade(name="E")
    new_root = Clade(name="root", clades=[tree, E])  # add E=A as outgroup
    sequences2 = {"A":"G", "B":"G", "C":"A", "D":"A", "E":"A"}
    rep2, sets2 = reconstruct_ancestral_sequence_with_sets(new_root, sequences2, new_root)
    mrca_set2 = sets2[0]
    assert mrca_set2 == {"A"}
    ctrl = ["A","A"]
    conv = ["G","G"]
    agree = sum(1 for r in ctrl if r in mrca_set2) / len(ctrl)
    assert agree >= 1.0
    # Count since convergent residue not in MRCA set
    assert all(r not in mrca_set2 for r in conv)
    assert conv.count("G") >= 2


def test_x_is_treated_as_missing_in_downpass():
    """Ensure X/x at tips are treated as missing and do not appear in MRCA sets."""
    tree = create_simple_tree()
    sequences = {"A":"X", "B":"G", "C":"A", "D":"A"}
    rep, sets = reconstruct_ancestral_sequence_with_sets(tree, sequences, tree)
    mrca_set = sets[0]
    # X should not be in MRCA set
    assert "X" not in mrca_set


# ─────────────────────────────────────────────────────────────────────────────
# Test integrated ancestral outgroup reconstruction
# ─────────────────────────────────────────────────────────────────────────────

def test_get_ancestral_outgroup_basic(outgroup_tree_file):
    """Test basic ancestral outgroup reconstruction."""
    # Tree: (((A,B)AB,(C,D)CD)ABCD,E)root
    # Analysis species: A, B, C, D
    # Alignment: all species present
    
    alignment_sequences = {
        "A": "GAT",
        "B": "GAT",
        "C": "GAT",
        "D": "GAT",
        "E": "CTG",
    }
    
    analysis_species = {"A", "B", "C", "D"}
    
    ancestral_seq, ancestral_id = get_ancestral_outgroup_for_alignment(
        outgroup_tree_file, alignment_sequences, analysis_species
    )
    
    # Should successfully reconstruct
    assert len(ancestral_seq) == 3
    assert ancestral_id == "ANCESTRAL_MRCA"


def test_get_ancestral_outgroup_subset_alignment(outgroup_tree_file):
    """Test reconstruction when alignment has subset of species."""
    # Tree: (((A,B)AB,(C,D)CD)ABCD,E)root
    # Analysis species: A, B, C, D
    # Alignment: only A, C, E present
    
    alignment_sequences = {
        "A": "GAT",
        "C": "GAC",
        "E": "CTG",
    }
    
    analysis_species = {"A", "B", "C", "D"}
    
    ancestral_seq, ancestral_id = get_ancestral_outgroup_for_alignment(
        outgroup_tree_file, alignment_sequences, analysis_species
    )
    
    # Should reconstruct for pruned tree containing A, C, E
    assert len(ancestral_seq) == 3
    assert ancestral_id == "ANCESTRAL_MRCA"


def test_get_ancestral_outgroup_mrca_at_root():
    """Test that reconstruction fails when MRCA is at root (no outgroup)."""
    # Tree without outgroup: ((A,B)AB,(C,D)CD)root
    tree = create_simple_tree()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        write_newick_to_file(tree, f.name)
        tree_file = f.name
    
    try:
        alignment_sequences = {
            "A": "GAT",
            "B": "GAT",
            "C": "GAT",
            "D": "GAT",
        }
        
        analysis_species = {"A", "B", "C", "D"}
        
        # Should raise error because MRCA is at root
        with pytest.raises(AncestralReconstructionError, match="MRCA.*at root"):
            get_ancestral_outgroup_for_alignment(
                tree_file, alignment_sequences, analysis_species
            )
    finally:
        os.unlink(tree_file)


def test_get_ancestral_outgroup_alignment_subset_mrca_at_root(outgroup_tree_file):
    """Test that reconstruction fails when pruned tree has MRCA at root."""
    # Tree: (((A,B)AB,(C,D)CD)ABCD,E)root
    # Analysis species: A, B, C, D
    # Alignment: only A, B, C, D (no E) -> pruned tree has no outgroup
    
    alignment_sequences = {
        "A": "GAT",
        "B": "GAT",
        "C": "GAT",
        "D": "GAT",
    }
    
    analysis_species = {"A", "B", "C", "D"}
    
    # Should raise error because E is missing, so MRCA is at pruned root
    with pytest.raises(AncestralReconstructionError, match="MRCA.*at root"):
        get_ancestral_outgroup_for_alignment(
            outgroup_tree_file, alignment_sequences, analysis_species
        )


def test_get_ancestral_outgroup_no_analysis_species_in_alignment(outgroup_tree_file):
    """Test error when no analysis species are in alignment."""
    alignment_sequences = {
        "E": "GAT",
    }
    
    analysis_species = {"A", "B", "C", "D"}
    
    with pytest.raises(AncestralReconstructionError, match="None of the analysis species"):
        get_ancestral_outgroup_for_alignment(
            outgroup_tree_file, alignment_sequences, analysis_species
        )


def test_get_ancestral_outgroup_partial_analysis_species(outgroup_tree_file):
    """Test reconstruction when only some analysis species are in alignment."""
    # Tree: (((A,B)AB,(C,D)CD)ABCD,E)root
    # Analysis species: A, B, C, D
    # Alignment: only A, B, E
    
    alignment_sequences = {
        "A": "GAT",
        "B": "GTC",
        "E": "CTG",
    }
    
    analysis_species = {"A", "B", "C", "D"}
    
    # Should successfully reconstruct using A and B (MRCA is AB clade)
    ancestral_seq, ancestral_id = get_ancestral_outgroup_for_alignment(
        outgroup_tree_file, alignment_sequences, analysis_species
    )
    
    assert len(ancestral_seq) == 3
    assert ancestral_id == "ANCESTRAL_MRCA"


# ─────────────────────────────────────────────────────────────────────────────
# Test tree validation for Fast Scan
# ─────────────────────────────────────────────────────────────────────────────

def test_validate_tree_valid(outgroup_tree_file, species_groups_file):
    """Test validation of a valid tree."""
    valid, message = validate_tree_for_fast_scan(outgroup_tree_file, species_groups_file)
    
    assert valid is True
    assert "valid" in message.lower()


def test_validate_tree_missing_species(outgroup_tree_file, tmp_path):
    """Test validation fails when tree missing analysis species."""
    # Create groups file with species not in tree
    groups_file = tmp_path / "bad_groups.txt"
    with open(groups_file, 'w') as f:
        f.write("X\n")
        f.write("Y\n")
        f.write("Z\n")
        f.write("W\n")
    
    valid, message = validate_tree_for_fast_scan(outgroup_tree_file, str(groups_file))
    
    assert valid is False
    assert "missing" in message.lower()


def test_validate_tree_mrca_at_root(simple_tree_file, species_groups_file):
    """Test validation fails when MRCA would be at root."""
    # simple_tree has no outgroup, so MRCA of all analysis species is root
    valid, message = validate_tree_for_fast_scan(simple_tree_file, species_groups_file)
    
    assert valid is False
    assert "root" in message.lower() or "outgroup" in message.lower()


def test_validate_tree_nonexistent_file(species_groups_file):
    """Test validation fails for nonexistent tree."""
    valid, message = validate_tree_for_fast_scan("/nonexistent.nwk", species_groups_file)
    
    assert valid is False
    assert "error" in message.lower() or "not found" in message.lower()


# ─────────────────────────────────────────────────────────────────────────────
# Test edge cases
# ─────────────────────────────────────────────────────────────────────────────

def test_prune_tree_preserves_original():
    """Test that pruning doesn't modify the original tree."""
    tree = create_simple_tree()
    original_terminals = get_terminal_names(tree)
    
    keep = {"A", "B"}
    pruned = prune_tree(tree, keep)
    
    # Original should be unchanged
    assert get_terminal_names(tree) == original_terminals
    # Pruned should have subset
    assert get_terminal_names(pruned) == keep


def test_reconstruct_empty_sequences():
    """Test reconstruction with empty sequences raises error."""
    tree = create_simple_tree()
    sequences = {}
    
    with pytest.raises(AncestralReconstructionError, match="No sequences"):
        reconstruct_ancestral_sequence(tree, sequences, tree)


def test_reconstruct_all_gaps():
    """Test reconstruction handles all-gap positions."""
    tree = create_simple_tree()
    sequences = {
        "A": "---",
        "B": "---",
        "C": "---",
        "D": "---",
    }
    
    ancestral_seq = reconstruct_ancestral_sequence(tree, sequences, tree)
    
    # Should return something (likely all '?')
    assert len(ancestral_seq) == 3


def test_fitch_parsimony_downpass_empty_children():
    """Test downpass handles nodes with no valid children."""
    tree = create_simple_tree()
    sequences = {
        "A": "-",
        "B": "-",
        "C": "G",
        "D": "T",
    }
    
    state_sets = fitch_parsimony_downpass(tree, sequences, 0)
    
    # AB clade should have empty set (both children are gaps)
    # CD clade should have {G, T}
    # Root should union them
    assert isinstance(state_sets, dict)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
