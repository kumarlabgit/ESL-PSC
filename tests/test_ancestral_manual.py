#!/usr/bin/env python3
"""Manual tests for ancestral reconstruction without pytest dependency."""

import sys
import os
import tempfile
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

# Add parent to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from gui.core.ancestral_reconstruction import (
    parse_tree,
    get_terminal_names,
    prune_tree,
    find_mrca,
    is_mrca_at_root,
    reconstruct_ancestral_sequence,
    get_ancestral_outgroup_for_alignment,
    validate_tree_for_fast_scan,
    AncestralReconstructionError,
)


def create_test_tree_with_outgroup():
    """Create tree: (((A,B)AB,(C,D)CD)ABCD,E)root."""
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


def write_newick(tree, path):
    """Write tree to file."""
    temp_tree = Tree(root=tree)
    Phylo.write(temp_tree, path, 'newick')


def test_basic_parsing():
    """Test tree parsing."""
    print("Testing tree parsing...")
    tree = create_test_tree_with_outgroup()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        write_newick(tree, f.name)
        tree_file = f.name
    
    try:
        parsed = parse_tree(tree_file)
        terminals = get_terminal_names(parsed)
        assert "A" in terminals and "E" in terminals
        print("  ✓ Tree parsing works")
    finally:
        os.unlink(tree_file)


def test_tree_pruning():
    """Test tree pruning."""
    print("Testing tree pruning...")
    tree = create_test_tree_with_outgroup()
    
    # Prune to just A, B, E
    pruned = prune_tree(tree, {"A", "B", "E"})
    terminals = get_terminal_names(pruned)
    
    assert terminals == {"A", "B", "E"}
    assert "C" not in terminals
    print("  ✓ Tree pruning works")


def test_mrca_finding():
    """Test MRCA identification."""
    print("Testing MRCA identification...")
    tree = create_test_tree_with_outgroup()
    
    # MRCA of A,B should be AB clade
    mrca_ab = find_mrca(tree, {"A", "B"})
    desc = get_terminal_names(mrca_ab)
    assert "A" in desc and "B" in desc
    assert "C" not in desc
    
    # MRCA of A,B,C,D should be ABCD clade (not root)
    mrca_abcd = find_mrca(tree, {"A", "B", "C", "D"})
    assert not is_mrca_at_root(tree, mrca_abcd)
    
    print("  ✓ MRCA identification works")


def test_parsimony_reconstruction():
    """Test Fitch parsimony reconstruction."""
    print("Testing parsimony reconstruction...")
    tree = create_test_tree_with_outgroup()
    
    sequences = {
        "A": "GATTACA",
        "B": "GATTACA",
        "C": "GATTAGA",
        "D": "GATTAGA",
        "E": "CTGCTGC",
    }
    
    # Reconstruct at ABCD node
    mrca_abcd = find_mrca(tree, {"A", "B", "C", "D"})
    ancestral = reconstruct_ancestral_sequence(tree, sequences, mrca_abcd)
    
    assert len(ancestral) == 7
    # First 5 positions identical across A,B,C,D
    assert ancestral[:5] == "GATTA"
    
    print("  ✓ Parsimony reconstruction works")


def test_integrated_reconstruction():
    """Test full integrated ancestral outgroup reconstruction."""
    print("Testing integrated reconstruction...")
    tree = create_test_tree_with_outgroup()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        write_newick(tree, f.name)
        tree_file = f.name
    
    try:
        alignment_sequences = {
            "A": "GAT",
            "B": "GAT",
            "C": "GAC",
            "E": "CTG",
        }
        
        analysis_species = {"A", "B", "C", "D"}
        
        ancestral_seq, ancestral_id = get_ancestral_outgroup_for_alignment(
            tree_file, alignment_sequences, analysis_species
        )
        
        assert len(ancestral_seq) == 3
        assert ancestral_id == "ANCESTRAL_MRCA"
        print("  ✓ Integrated reconstruction works")
    finally:
        os.unlink(tree_file)


def test_mrca_at_root_detection():
    """Test that MRCA at root is properly detected and rejected."""
    print("Testing MRCA at root detection...")
    
    # Create tree without outgroup
    A = Clade(name="A")
    B = Clade(name="B")
    C = Clade(name="C")
    D = Clade(name="D")
    AB = Clade(name="AB", clades=[A, B])
    CD = Clade(name="CD", clades=[C, D])
    root = Clade(name="root", clades=[AB, CD])
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        write_newick(root, f.name)
        tree_file = f.name
    
    try:
        alignment_sequences = {
            "A": "GAT",
            "B": "GAT",
            "C": "GAT",
            "D": "GAT",
        }
        
        analysis_species = {"A", "B", "C", "D"}
        
        # Should raise error
        try:
            get_ancestral_outgroup_for_alignment(
                tree_file, alignment_sequences, analysis_species
            )
            assert False, "Should have raised error"
        except AncestralReconstructionError as e:
            assert "root" in str(e).lower()
            print("  ✓ MRCA at root correctly detected")
    finally:
        os.unlink(tree_file)


def test_validation():
    """Test tree validation."""
    print("Testing tree validation...")
    tree = create_test_tree_with_outgroup()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        write_newick(tree, f.name)
        tree_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("A\nC\nB\nD\n")
        groups_file = f.name
    
    try:
        valid, message = validate_tree_for_fast_scan(tree_file, groups_file)
        assert valid, f"Validation failed: {message}"
        print("  ✓ Tree validation works")
    finally:
        os.unlink(tree_file)
        os.unlink(groups_file)


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("Running Ancestral Reconstruction Tests")
    print("="*60 + "\n")
    
    tests = [
        test_basic_parsing,
        test_tree_pruning,
        test_mrca_finding,
        test_parsimony_reconstruction,
        test_integrated_reconstruction,
        test_mrca_at_root_detection,
        test_validation,
    ]
    
    failed = 0
    for test in tests:
        try:
            test()
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    if failed == 0:
        print("✓ All tests passed!")
    else:
        print(f"✗ {failed} test(s) failed")
    print("="*60 + "\n")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
