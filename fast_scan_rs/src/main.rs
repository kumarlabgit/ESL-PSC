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

            let out_seq = species_seq.get(&outgroup).cloned().unwrap_or_else(|| String::new());

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

                    let clean_conv: Vec<char> = conv_ns.iter().copied().filter(|c| *c != '?' && *c != '-').collect();
                    let clean_ctrl: Vec<char> = ctrl_ns.iter().copied().filter(|c| *c != '?' && *c != '-').collect();
                    let clean_out: Vec<char> = out_aa.iter().copied().filter(|c| *c != '?' && *c != '-').collect();

                    if eligible_true {
                        if !clean_ctrl.is_empty() && !clean_out.is_empty() {
                            if let (Some(c_ctrl), Some(c_out)) = (is_uniform_non_gap(&clean_ctrl), is_uniform_non_gap(&clean_out)) {
                                if c_ctrl == c_out {
                                    let mut cnt: HashMap<char, u32> = HashMap::new();
                                    for r in clean_conv.iter() { *cnt.entry(*r).or_insert(0) += 1; }
                                    if cnt.iter().any(|(r, c)| *r != c_ctrl && *c >= 2) {
                                        ccs += 1;
                                    }
                                }
                            }
                        }
                    }

                    if eligible_ctrl {
                        if let Some(out_res) = is_uniform_non_gap(&clean_out) {
                            if !clean_conv.is_empty() && clean_conv.iter().all(|r| *r == out_res) {
                                let mut cnt: HashMap<char, u32> = HashMap::new();
                                for r in clean_ctrl.iter() { *cnt.entry(*r).or_insert(0) += 1; }
                                if cnt.iter().any(|(r, c)| *r != out_res && *c >= 2) {
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
