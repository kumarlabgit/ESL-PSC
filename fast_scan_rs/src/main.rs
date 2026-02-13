use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::io::{self, Read, Write};
use std::path::Path;
use rayon::prelude::*;
use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};

#[derive(Deserialize)]
struct Combo {
    conv: Vec<String>,
    ctrl: Vec<String>,
}

#[derive(Deserialize)]
struct InputSpec {
    alignment_dir: String,
    files: Option<Vec<String>>, // if None, scan all .fas
    combos: Vec<Combo>,
    outgroup: String,
    cs_threshold: Option<u32>,
    emit_progress: Option<bool>,
    min_out_ctrl_agreement: Option<f64>,
    tree_file: Option<String>,  // Newick tree for ancestral reconstruction
    analysis_species: Option<Vec<String>>,  // Species in analysis (for MRCA)
    tree_json: Option<NodeJson>, // Parsed tree provided by Python (preferred)
    require_unambiguous_mrca: Option<bool>, // Require single-residue MRCA (default false)
    compute_mrca_representative: Option<bool>, // Compute representative MRCA sequence (default false)
}

#[derive(Serialize)]
struct ResultRow {
    gene: String,
    avg_true: f64,
    avg_control: f64,
    diff: f64,
    variable_sites: u32,
    cs_sites_ge_4: u32,
    k_pairs: u32,
    per_combo_true: Vec<Option<f64>>,
    per_combo_diff: Vec<Option<f64>>,
}

fn parse_fasta(path: &Path) -> HashMap<String, String> {
    let mut map: HashMap<String, String> = HashMap::new();
    let Ok(text) = fs::read_to_string(path) else { return map };
    let mut cur_id: Option<String> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            cur_id = Some(rest.trim().to_string());
            map.entry(cur_id.clone().unwrap()).or_insert_with(String::new);
        } else {
            if let Some(id) = &cur_id {
                map.entry(id.clone()).and_modify(|s| s.push_str(line.trim()));
            }
        }
    }
    map
}

fn get_char(seq: &str, pos: usize) -> char {
    if pos >= seq.len() { return '?'; }
    // sequences are ASCII AAs; index by bytes is safe for ASCII
    seq.as_bytes()[pos] as char
}

fn is_uniform_non_gap(list: &[char]) -> Option<char> {
    let mut it = list.iter().copied().filter(|c| *c != '?' && *c != '-');
    if let Some(first) = it.next() {
        if it.clone().all(|c| c == first) {
            return Some(first);
        }
    }
    None
}

// ===== Tree and Ancestral Reconstruction =====

#[derive(Clone, Deserialize)]
struct NodeJson {
    name: Option<String>,
    children: Vec<NodeJson>,
}

#[derive(Clone)]
struct TreeNode {
    name: Option<String>,
    children: Vec<TreeNode>,
}

impl TreeNode {
    fn new(name: Option<String>) -> Self {
        TreeNode { name, children: vec![] }
    }
    
    fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
    
    fn get_terminals(&self) -> Vec<String> {
        if self.is_leaf() {
            if let Some(name) = &self.name {
                vec![name.clone()]
            } else {
                vec![]
            }
        } else {
            self.children.iter().flat_map(|c| c.get_terminals()).collect()
        }
    }
}

fn tree_from_json(n: &NodeJson) -> TreeNode {
    TreeNode {
        name: n.name.clone(),
        children: n.children.iter().map(tree_from_json).collect(),
    }
}

fn parse_newick(text: &str) -> Option<TreeNode> {
    let text = text.trim().trim_end_matches(';');
    parse_newick_node(text).map(|(node, _)| node)
}

fn parse_newick_node(s: &str) -> Option<(TreeNode, &str)> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    
    if s.starts_with('(') {
        // Internal node with children
        let s = &s[1..]; // Skip '('
        let mut children = vec![];
        let mut remaining = s;
        
        loop {
            if let Some((child, rest)) = parse_newick_node(remaining) {
                children.push(child);
                remaining = rest.trim();
                
                if remaining.starts_with(',') {
                    remaining = &remaining[1..];
                } else if remaining.starts_with(')') {
                    remaining = &remaining[1..];
                    break;
                } else {
                    return None;
                }
            } else {
                return None;
            }
        }
        
        // Parse node label (optional) and branch length (optional)
        let (name, remaining) = parse_node_label_and_length(remaining);
        let mut node = TreeNode::new(name);
        node.children = children;
        Some((node, remaining))
    } else {
        // Leaf node
        let (name, remaining) = parse_node_label_and_length(s);
        Some((TreeNode::new(name), remaining))
    }
}

fn parse_node_label_and_length(s: &str) -> (Option<String>, &str) {
    let s = s.trim();
    
    // Find where label ends (at ',' or ')' or ':' or end)
    let mut label_end = 0;
    let mut in_name = true;
    let chars: Vec<char> = s.chars().collect();
    
    for (i, &c) in chars.iter().enumerate() {
        if c == ':' || c == ',' || c == ')' {
            label_end = i;
            in_name = false;
            break;
        }
    }
    
    if in_name {
        label_end = chars.len();
    }
    
    let label = if label_end > 0 {
        Some(s[..label_end].trim().to_string())
    } else {
        None
    };
    
    // Skip branch length if present (after ':')
    let mut remaining = &s[label_end..];
    if remaining.starts_with(':') {
        // Find end of branch length
        if let Some(pos) = remaining[1..].find(|c| c == ',' || c == ')') {
            remaining = &remaining[pos + 1..];
        } else {
            remaining = "";
        }
    }
    
    (label, remaining)
}

fn prune_tree(node: &TreeNode, keep_species: &[String]) -> TreeNode {
    if node.is_leaf() {
        // Keep leaf if in keep_species
        if let Some(name) = &node.name {
            if keep_species.contains(name) {
                return node.clone();
            }
        }
        // Return empty node (will be removed)
        TreeNode::new(None)
    } else {
        // Recursively prune children
        let mut pruned_children: Vec<TreeNode> = node.children
            .iter()
            .map(|c| prune_tree(c, keep_species))
            .filter(|c| !c.children.is_empty() || c.name.is_some())
            .collect();
        
        // Collapse single-child nodes
        while pruned_children.len() == 1 && !pruned_children[0].is_leaf() {
            let child = pruned_children.remove(0);
            pruned_children = child.children;
        }
        
        let mut result = TreeNode::new(node.name.clone());
        result.children = pruned_children;
        result
    }
}

fn find_mrca<'a>(node: &'a TreeNode, species: &[String]) -> Option<&'a TreeNode> {
    let terminals = node.get_terminals();
    let has_all = species.iter().all(|s| terminals.contains(s));
    
    if !has_all {
        return None;
    }
    
    // Check if any child contains all species
    for child in &node.children {
        if let Some(mrca) = find_mrca(child, species) {
            return Some(mrca);
        }
    }
    
    // This node is the MRCA
    Some(node)
}

/// Reconstruct ancestral sequences, returning both a representative sequence
/// and the full state sets for ambiguous positions.
fn reconstruct_ancestral_with_sets(
    node: &TreeNode,
    sequences: &HashMap<String, String>,
    seq_len: usize
) -> (String, Vec<Vec<char>>) {
    let mut ancestral = String::new();
    let mut state_sets = Vec::new();
    
    for pos in 0..seq_len {
        // Downpass: compute state sets for each node
        let state_set = downpass(node, sequences, pos);
        
        // Pick representative (lexicographically first) for the ancestral sequence string
        let anc_char = if !state_set.is_empty() {
            let mut chars: Vec<char> = state_set.clone().into_iter().collect();
            chars.sort();
            chars[0]
        } else {
            '?'
        };
        
        ancestral.push(anc_char);
        state_sets.push(state_set);
    }
    
    (ancestral, state_sets)
}

fn node_key(node: &TreeNode) -> usize {
    node as *const TreeNode as usize
}

fn sort_dedup(mut v: Vec<char>) -> Vec<char> {
    if v.len() > 1 {
        v.sort();
        v.dedup();
    }
    v
}

fn downpass_collect(
    node: &TreeNode,
    sequences: &HashMap<String, String>,
    pos: usize,
    out_sets: &mut HashMap<usize, Vec<char>>,
) -> Vec<char> {
    let key = node_key(node);
    if node.is_leaf() {
        let set = if let Some(name) = &node.name {
            if let Some(seq) = sequences.get(name) {
                let ch = get_char(seq, pos);
                if ch != '-' && ch != '?' && ch != 'X' && ch != 'x' {
                    vec![ch]
                } else {
                    vec![]
                }
            } else {
                vec![]
            }
        } else {
            vec![]
        };
        out_sets.insert(key, set.clone());
        return set;
    }

    let mut child_sets: Vec<Vec<char>> = Vec::new();
    for c in &node.children {
        let s = downpass_collect(c, sequences, pos, out_sets);
        if !s.is_empty() {
            child_sets.push(s);
        }
    }
    if child_sets.is_empty() {
        out_sets.insert(key, vec![]);
        return vec![];
    }

    let mut result = child_sets[0].clone();
    for set in &child_sets[1..] {
        result.retain(|c| set.contains(c));
    }
    if result.is_empty() {
        // Union
        for set in &child_sets {
            for &c in set {
                if !result.contains(&c) {
                    result.push(c);
                }
            }
        }
    }
    let result = sort_dedup(result);
    out_sets.insert(key, result.clone());
    result
}

fn uppass_assign(
    node: &TreeNode,
    parent_state: Option<char>,
    sets: &HashMap<usize, Vec<char>>,
    out_assign: &mut HashMap<usize, char>,
) {
    let key = node_key(node);
    let node_set = sets.get(&key).cloned().unwrap_or_default();
    let state = if node_set.is_empty() {
        '?'
    } else if let Some(p) = parent_state {
        if node_set.contains(&p) {
            p
        } else {
            *node_set.iter().min().unwrap()
        }
    } else {
        *node_set.iter().min().unwrap()
    };
    out_assign.insert(key, state);
    for c in &node.children {
        uppass_assign(c, Some(state), sets, out_assign);
    }
}

/// Reconstruct the MRCA representative sequence using a Fitch uppass, while
/// also returning the MRCA downpass state sets per position.
///
/// This matches the Python implementation semantics:
/// - downpass computes state sets
/// - root assignment picks lexicographically-first state
/// - children inherit parent state if it is in their set, else pick lexicographically-first
fn reconstruct_mrca_with_sets(
    root: &TreeNode,
    mrca_key: usize,
    sequences: &HashMap<String, String>,
    seq_len: usize,
    compute_representative: bool,
) -> (String, Vec<Vec<char>>) {
    let mut ancestral = String::new();
    let mut mrca_state_sets: Vec<Vec<char>> = Vec::new();
    for pos in 0..seq_len {
        let mut sets: HashMap<usize, Vec<char>> = HashMap::new();
        let _root_set = downpass_collect(root, sequences, pos, &mut sets);

        let anc_char = if compute_representative {
            let mut assign: HashMap<usize, char> = HashMap::new();
            uppass_assign(root, None, &sets, &mut assign);
            assign.get(&mrca_key).copied().unwrap_or('?')
        } else {
            '?'
        };
        ancestral.push(anc_char);
        mrca_state_sets.push(sets.get(&mrca_key).cloned().unwrap_or_default());
    }
    (ancestral, mrca_state_sets)
}

/// Legacy wrapper for backward compatibility
fn reconstruct_ancestral(node: &TreeNode, sequences: &HashMap<String, String>, seq_len: usize) -> String {
    reconstruct_ancestral_with_sets(node, sequences, seq_len).0
}

fn downpass(node: &TreeNode, sequences: &HashMap<String, String>, pos: usize) -> Vec<char> {
    if node.is_leaf() {
        // Leaf: return observed state (treat X/x as missing like - and ?)
        if let Some(name) = &node.name {
            if let Some(seq) = sequences.get(name) {
                let ch = get_char(seq, pos);
                if ch != '-' && ch != '?' && ch != 'X' && ch != 'x' {
                    return vec![ch];
                }
            }
        }
        return vec![];
    }
    
    // Internal node: combine children
    let child_sets: Vec<Vec<char>> = node.children
        .iter()
        .map(|c| downpass(c, sequences, pos))
        .filter(|s| !s.is_empty())
        .collect();
    
    if child_sets.is_empty() {
        return vec![];
    }
    
    // Find intersection
    let mut result: Vec<char> = child_sets[0].clone();
    for set in &child_sets[1..] {
        result.retain(|c| set.contains(c));
    }
    
    if result.is_empty() {
        // Union
        for set in &child_sets {
            for &c in set {
                if !result.contains(&c) {
                    result.push(c);
                }
            }
        }
    }
    
    result
}

fn main() {
    // Read JSON from stdin
    let mut buf = String::new();
    if io::stdin().read_to_string(&mut buf).is_err() {
        eprintln!("failed to read stdin");
        std::process::exit(1);
    }
    let spec: InputSpec = match serde_json::from_str(&buf) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("invalid json: {}", e);
            std::process::exit(2);
        }
    };

    let cs_threshold = spec.cs_threshold.unwrap_or(4);
    let dir = Path::new(&spec.alignment_dir);
    let files: Vec<String> = if let Some(v) = spec.files {
        v
    } else {
        let mut v: Vec<String> = Vec::new();
        if let Ok(rd) = fs::read_dir(dir) {
            for ent in rd.flatten() {
                let p = ent.path();
                if let Some(ext) = p.extension() {
                    if ext == "fas" { if let Some(name) = p.file_name().and_then(|s| s.to_str()) { v.push(name.to_string()); } }
                }
            }
        }
        v.sort();
        v
    };

    let emit_progress = spec.emit_progress.unwrap_or(false);
    let progress = Arc::new(AtomicUsize::new(0));
    let combos_arc = Arc::new(spec.combos);
    let outgroup = spec.outgroup.clone();
    let mut min_agree = spec.min_out_ctrl_agreement.unwrap_or(1.0);
    if min_agree < 0.0 { min_agree = 0.0; }
    if min_agree > 1.0 { min_agree = 1.0; }
    
    let require_unambiguous_mrca = spec.require_unambiguous_mrca.unwrap_or(false);
    let compute_mrca_representative = spec.compute_mrca_representative.unwrap_or(false);
    
    // Build the tree once for this run. Prefer JSON supplied by Python.
    let tree_root = if let Some(ref j) = spec.tree_json {
        Some(tree_from_json(j))
    } else if let Some(ref tree_path) = spec.tree_file {
        if let Ok(tree_text) = fs::read_to_string(tree_path) {
            parse_newick(&tree_text)
        } else {
            None
        }
    } else {
        None
    };
    let tree_arc = Arc::new(tree_root);
    let analysis_species_arc = Arc::new(spec.analysis_species.clone().unwrap_or_default());

    let total_files = files.len();
    let results: Vec<ResultRow> = files
        .par_iter()
        .map(|fname| {
            let path = dir.join(fname);
            let species_seq = parse_fasta(&path);
            if species_seq.is_empty() {
                if emit_progress {
                    let cur = progress.fetch_add(1, Ordering::SeqCst) + 1;
                    if cur % 200 == 0 || cur == total_files {
                        let _ = writeln!(io::stderr(), "PROGRESS {}", cur);
                        let _ = io::stderr().flush();
                    }
                }
                return ResultRow{
                    gene: fname.trim_end_matches(".fas").to_string(),
                    avg_true: 0.0,
                    avg_control: 0.0,
                    diff: 0.0,
                    variable_sites: 0,
                    cs_sites_ge_4: 0,
                    k_pairs: 0,
                    per_combo_true: vec![],
                    per_combo_diff: vec![],
                };
            }
            let seq_len = species_seq.values().map(|s| s.len()).max().unwrap_or(0);
            // variable sites across all species
            let all_species: Vec<&String> = species_seq.keys().collect();
            let mut variable_sites_total: u32 = 0;
            for pos in 0..seq_len {
                let mut counts: HashMap<char, u32> = HashMap::new();
                for sp in &all_species {
                    let aa = get_char(species_seq.get(*sp).unwrap(), pos);
                    *counts.entry(aa).or_insert(0) += 1;
                }
                counts.remove(&'-');
                let num_left: u32 = counts.values().sum();
                if counts.len() <= 1 { continue; }
                if (counts.len() as u32) == num_left { continue; }
                if counts.values().max().copied().unwrap_or(0) == num_left - 1 { continue; }
                variable_sites_total += 1;
            }

            let mut files_true_counts: Vec<u32> = Vec::new();
            let mut files_ctrl_counts: Vec<u32> = Vec::new();
            let mut true_den: u32 = 0;
            let mut ctrl_den: u32 = 0;

            let mut best_cs_sites: u32 = 0;
            let mut best_k_pairs: u32 = 0;

            let mut per_combo_true: Vec<Option<f64>> = vec![None; combos_arc.len()];
            let mut per_combo_diff: Vec<Option<f64>> = vec![None; combos_arc.len()];

            // Compute outgroup sequence and MRCA state sets (for ancestral reconstruction)
            let (out_seq, mrca_state_sets): (String, Vec<Vec<char>>) = if let Some(ref tree) = *tree_arc {
                // Ancestral reconstruction
                if !analysis_species_arc.is_empty() {
                    // Get species present in this alignment
                    let alignment_species: Vec<String> = species_seq.keys()
                        .filter(|s| analysis_species_arc.contains(s))
                        .cloned()
                        .collect();
                    
                    if alignment_species.is_empty() {
                        (String::new(), vec![])  // Skip if no analysis species
                    } else {
                        // Prune tree to alignment species
                        let mut all_species_in_alignment: Vec<String> = species_seq.keys().cloned().collect();
                        all_species_in_alignment.sort();
                        let pruned = prune_tree(tree, &all_species_in_alignment);
                        
                        // Find MRCA of analysis species
                        if let Some(mrca) = find_mrca(&pruned, &alignment_species) {
                            // Check if MRCA is at root (would mean no outgroup)
                            let mrca_terminals = mrca.get_terminals();
                            let pruned_terminals = pruned.get_terminals();
                            if mrca_terminals.len() == pruned_terminals.len() {
                                (String::new(), vec![])  // Skip - MRCA at root
                            } else {
                                // Reconstruct MRCA representative using uppass on the full pruned tree
                                // (matches Python), while keeping downpass state sets for MRCA.
                                let mk = node_key(mrca);
                                reconstruct_mrca_with_sets(
                                    &pruned,
                                    mk,
                                    &species_seq,
                                    seq_len,
                                    compute_mrca_representative,
                                )
                            }
                        } else {
                            (String::new(), vec![])  // Skip if MRCA not found
                        }
                    }
                } else {
                    let seq = species_seq.get(&outgroup).cloned().unwrap_or_else(|| String::new());
                    (seq, vec![])
                }
            } else {
                // Use provided outgroup species
                let seq = species_seq.get(&outgroup).cloned().unwrap_or_else(|| String::new());
                (seq, vec![])
            };
            
            let use_mrca_sets = !mrca_state_sets.is_empty();
            let require_unamb = require_unambiguous_mrca;

            for (combo_idx, combo) in combos_arc.iter().enumerate() {
                let pairs_len = std::cmp::min(combo.conv.len(), combo.ctrl.len());
                let mut pairs_present: Vec<(&str, &str)> = Vec::new();
                for i in 0..pairs_len {
                    let a = combo.conv[i].as_str();
                    let b = combo.ctrl[i].as_str();
                    if species_seq.contains_key(a) && species_seq.contains_key(b) {
                        pairs_present.push((a, b));
                    }
                }
                let conv_present: Vec<&str> = pairs_present.iter().map(|(a, _)| *a).collect();
                let ctrl_present: Vec<&str> = pairs_present.iter().map(|(_, b)| *b).collect();
                let eligible_true = conv_present.len() >= 2 && ctrl_present.len() >= 1;
                let eligible_ctrl = ctrl_present.len() >= 2 && conv_present.len() >= 1;
                if !(eligible_true || eligible_ctrl) { continue; }

                let mut ccs: u32 = 0;
                let mut ctrl_conv: u32 = 0;
                let mut cs_sites: u32 = 0;

                for pos in 0..seq_len {
                    let conv_aa: Vec<char> = conv_present
                        .iter()
                        .map(|sp| get_char(species_seq.get(*sp).unwrap(), pos))
                        .collect();
                    let ctrl_aa: Vec<char> = ctrl_present
                        .iter()
                        .map(|sp| get_char(species_seq.get(*sp).unwrap(), pos))
                        .collect();
                    let out_aa: Vec<char> = if pos < out_seq.len() { vec![get_char(&out_seq, pos)] } else { vec!['?'] };

                    let mut cc_counts: HashMap<char, u32> = HashMap::new();
                    for c in conv_aa.iter().chain(ctrl_aa.iter()) {
                        *cc_counts.entry(*c).or_insert(0) += 1;
                    }
                    let mut conv_ns = conv_aa.clone();
                    let mut ctrl_ns = ctrl_aa.clone();
                    for lst in [&mut conv_ns, &mut ctrl_ns] {
                        for i in 0..lst.len() {
                            let r = lst[i];
                            if r != '-' && *cc_counts.get(&r).unwrap_or(&0) == 1 { lst[i] = '?'; }
                        }
                    }

                    let clean_conv: Vec<char> = conv_ns.iter().copied().filter(|c| *c != '?' && *c != '-' && *c != 'X' && *c != 'x').collect();
                    let clean_ctrl: Vec<char> = ctrl_ns.iter().copied().filter(|c| *c != '?' && *c != '-' && *c != 'X' && *c != 'x').collect();
                    let clean_out: Vec<char> = out_aa.iter().copied().filter(|c| *c != '?' && *c != '-' && *c != 'X' && *c != 'x').collect();

                    if eligible_true {
                        // Get MRCA possible states for this position (if using ancestral reconstruction)
                        let mrca_set: Vec<char> = if use_mrca_sets && pos < mrca_state_sets.len() {
                            mrca_state_sets[pos].clone()
                        } else if !clean_out.is_empty() {
                            vec![clean_out[0]]  // Single outgroup residue
                        } else {
                            vec![]  // No valid outgroup
                        };
                        
                        // Skip if MRCA is ambiguous and user requires unambiguous
                        if !mrca_set.is_empty() && !(require_unamb && mrca_set.len() > 1) {
                            if !clean_ctrl.is_empty() {
                                // Check if controls match any MRCA possible state
                                let ctrl_mrca_matches: Vec<char> = clean_ctrl.iter()
                                    .copied()
                                    .filter(|c| mrca_set.contains(c))
                                    .collect();
                                let matches: f64 = ctrl_mrca_matches.len() as f64;
                                let total: f64 = clean_ctrl.len() as f64;
                                let agree: f64 = if total > 0.0 { matches / total } else { 0.0 };
                                
                                if agree + f64::EPSILON >= min_agree {
                                    // Count CCS: convergent share a residue NOT in MRCA set
                                    let mut cnt: HashMap<char, u32> = HashMap::new();
                                    for r in clean_conv.iter() { *cnt.entry(*r).or_insert(0) += 1; }
                                    if cnt.iter().any(|(r, c)| !mrca_set.contains(r) && *c >= 2) {
                                        ccs += 1;
                                    }
                                }
                            }
                        }
                    }

                    if eligible_ctrl {
                        // Get MRCA possible states for this position (if using ancestral reconstruction)
                        let mrca_set: Vec<char> = if use_mrca_sets && pos < mrca_state_sets.len() {
                            mrca_state_sets[pos].clone()
                        } else if !clean_out.is_empty() {
                            vec![clean_out[0]]  // Single outgroup residue
                        } else {
                            vec![]  // No valid outgroup
                        };
                        
                        // Skip if MRCA is ambiguous and user requires unambiguous
                        if !mrca_set.is_empty() && !(require_unamb && mrca_set.len() > 1) {
                            // Check if all convergent match any MRCA possible state
                            if !clean_conv.is_empty() && clean_conv.iter().all(|r| mrca_set.contains(r)) {
                                let mut cnt: HashMap<char, u32> = HashMap::new();
                                for r in clean_ctrl.iter() { *cnt.entry(*r).or_insert(0) += 1; }
                                // Control convergence: control share a residue NOT in MRCA set
                                if cnt.iter().any(|(r, c)| !mrca_set.contains(r) && *c >= 2) {
                                    ctrl_conv += 1;
                                }
                            }
                        }
                    }

                    let mut diff_map: HashMap<char, i32> = HashMap::new();
                    for r in clean_conv.iter() { *diff_map.entry(*r).or_insert(0) += 1; }
                    for r in clean_ctrl.iter() { *diff_map.entry(*r).or_insert(0) -= 1; }
                    let raw_score: i32 = diff_map.values().map(|v| v.abs()).sum();
                    let gap_count: i32 = conv_ns.iter().filter(|c| **c == '-').count() as i32 + ctrl_ns.iter().filter(|c| **c == '-').count() as i32;
                    let mut final_score = raw_score - gap_count;
                    if final_score < 0 { final_score = 0; }
                    if (final_score as u32) >= cs_threshold {
                        cs_sites += 1;
                    }
                }

                if eligible_true {
                    files_true_counts.push(ccs);
                    true_den += 1;
                    per_combo_true[combo_idx] = Some(ccs as f64);
                }
                if eligible_ctrl {
                    files_ctrl_counts.push(ctrl_conv);
                    ctrl_den += 1;
                }
                if eligible_true && eligible_ctrl {
                    per_combo_diff[combo_idx] = Some((ccs as i32 - ctrl_conv as i32) as f64);
                }
                if cs_sites > best_cs_sites {
                    best_cs_sites = cs_sites;
                    best_k_pairs = pairs_present.len() as u32;
                }
            }

            let avg_true: f64 = if true_den > 0 { (files_true_counts.iter().sum::<u32>() as f64) / (true_den as f64) } else { 0.0 };
            let avg_ctrl: f64 = if ctrl_den > 0 { (files_ctrl_counts.iter().sum::<u32>() as f64) / (ctrl_den as f64) } else { 0.0 };

            let row = ResultRow{
                gene: fname.trim_end_matches(".fas").to_string(),
                avg_true,
                avg_control: avg_ctrl,
                diff: avg_true - avg_ctrl,
                variable_sites: variable_sites_total,
                cs_sites_ge_4: best_cs_sites,
                k_pairs: best_k_pairs,
                per_combo_true,
                per_combo_diff,
            };

            if emit_progress {
                let cur = progress.fetch_add(1, Ordering::SeqCst) + 1;
                if cur % 200 == 0 || cur == total_files {
                    let _ = writeln!(io::stderr(), "PROGRESS {}", cur);
                    let _ = io::stderr().flush();
                }
            }
            row
        })
        .collect();

    match serde_json::to_string(&results) {
        Ok(s) => {
            println!("{}", s);
        }
        Err(e) => {
            eprintln!("serialization error: {}", e);
            std::process::exit(3);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mrca_assignment_prefers_parent_state_like_python() {
        // root has children: mrca_subtree (A,B) and outgroup E.
        // At pos 0:
        // - MRCA downpass set is {A,T} (A leaf=A, B leaf=T)
        // - outgroup is T
        // => root set intersection is {T}, so root assignment is T
        // Python uppass would assign MRCA = T (inherit parent), not A (lexicographic min at MRCA).
        let mrca = TreeNode {
            name: Some("AB".to_string()),
            children: vec![
                TreeNode::new(Some("A".to_string())),
                TreeNode::new(Some("B".to_string())),
            ],
        };
        let root = TreeNode {
            name: Some("root".to_string()),
            children: vec![
                mrca.clone(),
                TreeNode::new(Some("E".to_string())),
            ],
        };

        let mut seqs: HashMap<String, String> = HashMap::new();
        seqs.insert("A".to_string(), "A".to_string());
        seqs.insert("B".to_string(), "T".to_string());
        seqs.insert("E".to_string(), "T".to_string());

        // Find MRCA node reference in this constructed tree (first child).
        let mrca_ref = &root.children[0];
        let mk = node_key(mrca_ref);
        let (rep, sets) = reconstruct_mrca_with_sets(&root, mk, &seqs, 1, true);
        assert_eq!(rep, "T");
        assert_eq!(sort_dedup(sets[0].clone()), vec!['A', 'T']);
    }
}
