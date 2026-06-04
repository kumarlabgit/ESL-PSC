"""Parsimony-based ancestral sequence reconstruction for Site Counter.

This module provides functionality to:
1. Parse and prune phylogenetic trees to match alignment species
2. Identify MRCA (Most Recent Common Ancestor) of analysis species
3. Perform Fitch parsimony reconstruction of ancestral sequences
4. Validate that reconstructions have appropriate outgroup context
"""
from __future__ import annotations

import os
from typing import Dict, List, Set, Tuple, Optional
from functools import lru_cache
from collections import Counter
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade


class AncestralReconstructionError(Exception):
    """Raised when ancestral reconstruction cannot proceed."""
    pass


def parse_tree(tree_path: str) -> Clade:
    """Parse a tree file (Newick or NEXUS).
    
    Parameters
    ----------
    tree_path : str
        Path to tree file
        
    Returns
    -------
    Clade
        Root clade of the parsed tree
        
    Raises
    ------
    AncestralReconstructionError
        If tree cannot be parsed
    """
    if not os.path.exists(tree_path):
        raise AncestralReconstructionError(f"Tree file not found: {tree_path}")
    
    try:
        def _looks_like_single_label_newick(s: str) -> bool:
            # Allow a single-leaf tree like "A;" (rare but technically Newick).
            # Reject arbitrary prose like "this is not a tree".
            import re
            s = s.strip()
            if s.endswith(";"):
                s = s[:-1].strip()
            return bool(re.fullmatch(r"[A-Za-z0-9_.-]+", s))

        # Try NEXUS first if extension suggests it
        ext = os.path.splitext(tree_path)[1].lower()
        if ext in {'.nexus', '.nex'}:
            tree = Phylo.read(tree_path, 'nexus')
            return tree.root
        
        # Check file content for NEXUS header
        with open(tree_path, 'r', encoding='utf-8', errors='ignore') as f:
            text = f.read()
            first_nonempty = ""
            for line in text.splitlines():
                if line.strip():
                    first_nonempty = line.strip()
                    break
            if first_nonempty.upper().startswith('#NEXUS'):
                tree = Phylo.read(tree_path, 'nexus')
                return tree.root

            # Quick sanity check for Newick-like content. Biopython will happily
            # parse arbitrary text as a one-leaf "tree" (leaf label), which is
            # not useful for Site Counter and breaks expectations in tests.
            stripped = text.strip()
            if not stripped:
                raise AncestralReconstructionError("Failed to parse tree file: empty file")
            # If there is no Newick punctuation, only allow a strict single-label tree.
            if not any(ch in stripped for ch in ("(", ")", ",", ":")):
                if not _looks_like_single_label_newick(stripped):
                    raise AncestralReconstructionError("Failed to parse tree file: not valid Newick or NEXUS")
        
        # Default to Newick
        tree = Phylo.read(tree_path, 'newick')
        return tree.root
    except Exception as e:
        raise AncestralReconstructionError(f"Failed to parse tree file: {e}")


def get_terminal_names(tree: Clade) -> Set[str]:
    """Get all terminal (leaf) names from a tree.
    
    Parameters
    ----------
    tree : Clade
        Root of tree or subtree
        
    Returns
    -------
    Set[str]
        Set of terminal names
    """
    names = set()
    for terminal in tree.get_terminals():
        if terminal.name:
            names.add(terminal.name)
    return names


def prune_tree(tree: Clade, keep_species: Set[str]) -> Clade:
    """Prune tree to only include specified species.
    
    Creates a deep copy of the tree and removes all terminals not in keep_species.
    Internal nodes with only one child after pruning are collapsed.
    
    Parameters
    ----------
    tree : Clade
        Root of the original tree
    keep_species : Set[str]
        Species to retain in the pruned tree
        
    Returns
    -------
    Clade
        Root of the pruned tree
        
    Raises
    ------
    AncestralReconstructionError
        If no species from keep_species are found in tree
    """
    import copy

    # Deep copy to avoid modifying original.
    # We then prune recursively so internal nodes that become terminal (e.g. "CD")
    # are dropped unless explicitly kept as a species label.
    root = copy.deepcopy(tree)

    def _prune_rec(node: Clade) -> Optional[Clade]:
        if node.is_terminal():
            return node if (node.name in keep_species) else None

        kept_children: List[Clade] = []
        for child in list(getattr(node, "clades", []) or []):
            kept = _prune_rec(child)
            if kept is not None:
                kept_children.append(kept)

        node.clades = kept_children

        # If an internal node ends up with no kept descendants, drop it (unless
        # its label is explicitly being used as a "species" identifier).
        if not node.clades:
            return node if (node.name in keep_species) else None

        # Collapse chains of single-child internal nodes to keep trees small and
        # avoid internal labels becoming terminals.
        if node.name not in keep_species:
            while (not node.is_terminal()) and len(node.clades) == 1 and node.name not in keep_species:
                child = node.clades[0]
                # Preserve branch length sums if present.
                if getattr(node, "branch_length", None) is not None:
                    if getattr(child, "branch_length", None) is None:
                        child.branch_length = node.branch_length
                    else:
                        child.branch_length += node.branch_length
                node = child

        return node

    pruned = _prune_rec(root)
    if pruned is None:
        terminals = get_terminal_names(tree)
        raise AncestralReconstructionError(
            f"None of the requested species found in tree. "
            f"Requested: {sorted(keep_species)[:5]}..., "
            f"Tree has: {sorted(terminals)[:5]}..."
        )

    return pruned


def _remove_empty_nodes(tree: Clade, keep_species: Set[str]) -> None:
    """Remove internal nodes that have no children (in-place).
    
    Only removes nodes that are not in the keep_species set (i.e., former
    internal nodes that became childless after pruning).
    
    Parameters
    ----------
    tree : Clade
        Root of tree to clean
    keep_species : Set[str]
        Species names to preserve
    """
    changed = True
    while changed:
        changed = False
        # Collect all nodes first to avoid modification during iteration
        all_nodes = list(tree.find_clades(order='postorder'))
        for node in all_nodes:
            # Skip root
            if node == tree:
                continue
            # If node has no children (became terminal after pruning) and is not
            # one of the species we want to keep, remove it
            if len(node.clades) == 0 and node.name not in keep_species:
                path = tree.get_path(node)
                if len(path) >= 2:
                    parent = path[-2]
                    parent.clades = [c for c in parent.clades if c != node]
                    changed = True
                    break


def _collapse_single_child_nodes(tree: Clade) -> None:
    """Collapse internal nodes that have only one child (in-place).
    
    Parameters
    ----------
    tree : Clade
        Root of tree to collapse
    """
    changed = True
    while changed:
        changed = False
        for node in list(tree.find_clades(order='level')):
            if node.is_terminal():
                continue
            # If node has exactly one child, collapse it
            if len(node.clades) == 1:
                child = node.clades[0]
                # Replace node's clades with child's clades
                node.clades = child.clades
                # Transfer any other attributes if needed
                if hasattr(child, 'branch_length') and child.branch_length is not None:
                    if hasattr(node, 'branch_length') and node.branch_length is not None:
                        node.branch_length += child.branch_length
                    else:
                        node.branch_length = child.branch_length
                changed = True
                break


def find_mrca(tree: Clade, species: Set[str]) -> Clade:
    """Find the Most Recent Common Ancestor of a set of species.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    species : Set[str]
        Set of species names
        
    Returns
    -------
    Clade
        The MRCA clade
        
    Raises
    ------
    AncestralReconstructionError
        If species not found or no common ancestor exists
    """
    if not species:
        raise AncestralReconstructionError("No species provided for MRCA search")
    
    # Get terminal clades matching species names
    terminals = []
    for terminal in tree.get_terminals():
        if terminal.name in species:
            terminals.append(terminal)
    
    if not terminals:
        raise AncestralReconstructionError(f"No matching species found in tree: {species}")
    
    if len(terminals) == 1:
        # Single species: return that terminal
        return terminals[0]
    
    # Use Bio.Phylo's common_ancestor method
    # Build a minimal Tree object for the method to work
    from Bio.Phylo.BaseTree import Tree
    temp_tree = Tree(root=tree)
    
    try:
        mrca = temp_tree.common_ancestor(terminals)
        return mrca
    except Exception as e:
        raise AncestralReconstructionError(f"Could not find MRCA: {e}")


def is_mrca_at_root(tree: Clade, mrca: Clade) -> bool:
    """Check if the MRCA is at the root of the tree.
    
    For Site Counter purposes, if the MRCA of analysis species is the root,
    there are no outgroup species in the alignment, so reconstruction is invalid.
    
    Parameters
    ----------
    tree : Clade
        Root of the tree
    mrca : Clade
        MRCA clade to check
        
    Returns
    -------
    bool
        True if MRCA is the root (after unwrapping single-child ancestors)
    """
    # A robust check for "no outgroup exists": if the MRCA spans all terminals
    # in the (possibly pruned) tree, then there are no terminals outside the
    # MRCA subtree.
    try:
        return get_terminal_names(tree) == get_terminal_names(mrca)
    except Exception:
        return tree == mrca


def fitch_parsimony_downpass(
    tree: Clade,
    sequences: Dict[str, str],
    position: int
) -> Dict[Clade, Set[str]]:
    """Perform Fitch parsimony downpass (tip to root) for a single alignment position.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    sequences : Dict[str, str]
        Map of species name to sequence
    position : int
        Alignment position (0-indexed)
        
    Returns
    -------
    Dict[Clade, Set[str]]
        Map of each clade to its possible residue set at this position
    """
    state_sets: Dict[Clade, Set[str]] = {}
    
    # Post-order traversal (children before parents)
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            # Leaf: use observed sequence
            if node.name and node.name in sequences:
                seq = sequences[node.name]
                if position < len(seq):
                    residue = seq[position]
                    # Skip gaps and missing data for parsimony
                    if residue not in {'-', '?', 'X', 'x'}:
                        state_sets[node] = {residue}
                    else:
                        state_sets[node] = set()
                else:
                    state_sets[node] = set()
            else:
                state_sets[node] = set()
        else:
            # Internal node: intersect or union child sets
            child_sets = [state_sets.get(child, set()) for child in node.clades]
            # Remove empty sets
            child_sets = [s for s in child_sets if s]
            
            if not child_sets:
                state_sets[node] = set()
            elif len(child_sets) == 1:
                state_sets[node] = child_sets[0].copy()
            else:
                # Fitch algorithm: intersection if non-empty, else union
                intersection = set.intersection(*child_sets)
                if intersection:
                    state_sets[node] = intersection
                else:
                    state_sets[node] = set.union(*child_sets)
    
    return state_sets


def fitch_parsimony_uppass(
    tree: Clade,
    state_sets: Dict[Clade, Set[str]]
) -> Dict[Clade, str]:
    """Perform Fitch parsimony uppass (root to tips) to assign final states.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    state_sets : Dict[Clade, Set[str]]
        State sets from downpass
        
    Returns
    -------
    Dict[Clade, str]
        Final assigned residue for each clade
    """
    assignments: Dict[Clade, str] = {}
    
    # Pre-order traversal (parents before children)
    for node in tree.find_clades(order='preorder'):
        node_set = state_sets.get(node, set())
        
        if not node_set:
            # No valid states: assign gap or missing
            assignments[node] = '?'
            continue
        
        if node == tree:
            # Root: pick most common state, or first alphabetically
            assignments[node] = min(node_set)
        else:
            # Non-root: prefer parent's assignment if in set, else pick first
            parent = _get_parent(tree, node)
            parent_state = assignments.get(parent, '?')
            if parent_state in node_set:
                assignments[node] = parent_state
            else:
                assignments[node] = min(node_set)
    
    return assignments


def _get_parent(tree: Clade, child: Clade) -> Optional[Clade]:
    """Get parent clade of a child in the tree.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    child : Clade
        Child clade
        
    Returns
    -------
    Optional[Clade]
        Parent clade, or None if child is root
    """
    path = tree.get_path(child)
    if len(path) >= 2:
        return path[-2]
    return None


def reconstruct_ancestral_sequence_with_sets(
    tree: Clade,
    sequences: Dict[str, str],
    mrca: Clade,
    *,
    compute_representative: bool = True,
) -> Tuple[str, List[Set[str]]]:
    """Reconstruct ancestral sequence at MRCA using Fitch parsimony, returning state sets.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    sequences : Dict[str, str]
        Map of species name to sequence
    mrca : Clade
        MRCA node to reconstruct
        
    Returns
    -------
    Tuple[str, List[Set[str]]]
        (ancestral_sequence, state_sets_per_position)
        ancestral_sequence: representative sequence (lexicographically first when ambiguous)
            when ``compute_representative=True``; otherwise placeholder '?' chars.
        state_sets_per_position: list of possible residue sets for each position
        
    Raises
    ------
    AncestralReconstructionError
        If reconstruction fails
    """
    if not sequences:
        raise AncestralReconstructionError("No sequences provided for reconstruction")
    
    # Determine alignment length
    seq_lengths = [len(s) for s in sequences.values()]
    if not seq_lengths:
        raise AncestralReconstructionError("All sequences are empty")
    
    max_length = max(seq_lengths)
    
    # Reconstruct position by position
    ancestral_seq = []
    mrca_state_sets = []
    for pos in range(max_length):
        # Downpass
        state_sets = fitch_parsimony_downpass(tree, sequences, pos)
        if compute_representative:
            # Uppass only needed when callers request a single-residue ancestral sequence.
            assignments = fitch_parsimony_uppass(tree, state_sets)
            residue = assignments.get(mrca, '?')
        else:
            residue = '?'
        mrca_set = state_sets.get(mrca, set())
        ancestral_seq.append(residue)
        mrca_state_sets.append(mrca_set)
    
    return ''.join(ancestral_seq), mrca_state_sets


def reconstruct_ancestral_sequence(
    tree: Clade,
    sequences: Dict[str, str],
    mrca: Clade
) -> str:
    """Reconstruct ancestral sequence at MRCA using Fitch parsimony.
    
    Legacy wrapper that returns only the representative sequence.
    
    Parameters
    ----------
    tree : Clade
        Root of tree
    sequences : Dict[str, str]
        Map of species name to sequence
    mrca : Clade
        MRCA node to reconstruct
        
    Returns
    -------
    str
        Reconstructed ancestral sequence
        
    Raises
    ------
    AncestralReconstructionError
        If reconstruction fails
    """
    return reconstruct_ancestral_sequence_with_sets(tree, sequences, mrca)[0]


@lru_cache(maxsize=8)
def _load_tree(tree_path: str) -> Clade:
    """Cached wrapper around ``parse_tree``.

    Keeps frequently-used trees resident in memory so Python fallback mode
    does not re-parse the same file for every alignment.
    """
    return parse_tree(tree_path)


def clade_to_json(node: Clade) -> dict:
    """Convert a Biopython Clade to a compact JSON-serializable dict.

    Schema:
      { "name": str | None, "children": [ ... ] }

    Parameters
    ----------
    node : Clade
        Root or subclade to serialize

    Returns
    -------
    dict
        JSON-serializable representation
    """
    return {
        "name": getattr(node, "name", None),
        "children": [clade_to_json(c) for c in getattr(node, "clades", [])],
    }


def get_ancestral_outgroup_for_alignment(
    tree_path: str,
    alignment_sequences: Dict[str, str],
    analysis_species: Set[str]
) -> Tuple[str, str]:
    """Get reconstructed ancestral outgroup sequence for an alignment.
    
    This is the main entry point for Site Counter integration.
    
    Parameters
    ----------
    tree_path : str
        Path to the full species tree (Newick or NEXUS)
    alignment_sequences : Dict[str, str]
        Sequences present in this alignment
    analysis_species : Set[str]
        Species involved in the convergence analysis (from species groups)
        
    Returns
    -------
    Tuple[str, str]
        (outgroup_sequence, identifier_for_logging)
        
    Raises
    ------
    AncestralReconstructionError
        If reconstruction cannot proceed (e.g., MRCA at root, no outgroup)
    """
    # Parse the tree (cached in Python fallback scenarios)
    full_tree = _load_tree(tree_path)
    
    # Identify MRCA of analysis species in the full tree
    analysis_in_tree = analysis_species & get_terminal_names(full_tree)
    if not analysis_in_tree:
        raise AncestralReconstructionError(
            f"None of the analysis species found in tree"
        )
    
    mrca_full = find_mrca(full_tree, analysis_in_tree)
    
    # Prune tree to alignment species
    alignment_species = set(alignment_sequences.keys())
    pruned_tree = prune_tree(full_tree, alignment_species)
    
    # Find the corresponding MRCA in the pruned tree
    # It should be the MRCA of the same analysis species that are present
    analysis_in_alignment = analysis_species & alignment_species
    if not analysis_in_alignment:
        raise AncestralReconstructionError(
            f"None of the analysis species found in this alignment"
        )
    
    mrca_pruned = find_mrca(pruned_tree, analysis_in_alignment)
    
    # Check if MRCA is at root of pruned tree
    if is_mrca_at_root(pruned_tree, mrca_pruned):
        raise AncestralReconstructionError(
            "MRCA of analysis species is at root of pruned tree; "
            "no outgroup species available for this alignment"
        )
    
    # Reconstruct ancestral sequence at MRCA
    ancestral_seq = reconstruct_ancestral_sequence(
        pruned_tree,
        alignment_sequences,
        mrca_pruned
    )
    
    # Generate a unique identifier for the reconstructed ancestor
    # Use a prefix that won't collide with real species names
    identifier = f"ANCESTRAL_MRCA"
    
    return ancestral_seq, identifier


def validate_tree_for_site_counter(
    tree_path: str,
    species_groups_file: str
) -> Tuple[bool, str]:
    """Validate that a tree is suitable for Site Counter ancestral reconstruction.
    
    Parameters
    ----------
    tree_path : str
        Path to tree file
    species_groups_file : str
        Path to species groups file
        
    Returns
    -------
    Tuple[bool, str]
        (is_valid, message)
    """
    try:
        # Parse tree
        tree = parse_tree(tree_path)
        
        # Basal bifurcation check (mirror InputPage.open_newick_tree)
        # Descend through any single-child wrappers from the root and require
        # that the first branching node has exactly two children. This avoids
        # relying on Bio.Phylo's rooted metadata.
        try:
            node = tree
            safety_counter = 0
            while (
                node is not None and hasattr(node, "clades") and isinstance(getattr(node, "clades", None), list)
                and len(node.clades) == 1 and safety_counter < 1000
            ):
                node = node.clades[0]
                safety_counter += 1
            base_children = len(node.clades) if node is not None and hasattr(node, "clades") else 0
            if base_children != 2:
                return False, (
                    "The tree does not appear to be properly rooted: the basal split is not a bifurcation. "
                    "Please root your tree before running Site Counter."
                )
        except Exception:
            # If detection fails, do not block; continue with other validations
            pass
        tree_species = get_terminal_names(tree)
        
        # Parse species groups
        from esl_psc_cli.deletion_canceler import parse_species_groups
        combos = parse_species_groups(species_groups_file)
        
        # Collect all analysis species
        analysis_species = set()
        for combo in combos:
            for sp in combo:
                analysis_species.add(sp)
        
        # Check coverage
        missing = analysis_species - tree_species
        if missing:
            return False, (
                f"Tree is missing {len(missing)} analysis species: "
                f"{sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}"
            )
        
        # Check that MRCA is not at root (would mean no outgroup)
        mrca = find_mrca(tree, analysis_species)
        if is_mrca_at_root(tree, mrca):
            return False, (
                "MRCA of all analysis species is at tree root; "
                "no outgroup available for ancestral reconstruction"
            )
        
        return True, "Tree is valid for ancestral reconstruction"
        
    except Exception as e:
        return False, f"Validation error: {e}"


# Backward compatibility alias.
validate_tree_for_fast_scan = validate_tree_for_site_counter
