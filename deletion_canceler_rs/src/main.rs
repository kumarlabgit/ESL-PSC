use clap::Parser;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser)]
#[command(author, version, about = "Rust deletion canceler", long_about = None)]
struct Args {
    #[arg(long)]
    alignments_dir: PathBuf,
    #[arg(long)]
    canceled_alignments_dir: Option<PathBuf>,
    #[arg(long, conflicts_with = "species_groups_file")]
    response_file: Option<PathBuf>,
    #[arg(long, conflicts_with = "response_file")]
    species_groups_file: Option<PathBuf>,
    #[arg(long, conflicts_with_all = &["species_groups_file", "response_file"])]
    response_dir: Option<PathBuf>,
    #[arg(long, default_value_t = false)]
    cancel_tri_allelic: bool,
    #[arg(long, default_value_t = false)]
    nix_full_deletions: bool,
    #[arg(long)]
    outgroup_species: Option<String>,
    #[arg(long, default_value_t = false)]
    cancel_only_partner: bool,
    #[arg(long, default_value_t = 2)]
    min_pairs: usize,
    #[arg(long)]
    limited_genes_list: Option<PathBuf>,
}

fn read_response(path: &Path) -> Result<Vec<String>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut species = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 && parts[1].trim() != "0" {
            species.push(parts[0].to_string());
        }
    }
    Ok(species)
}

fn parse_species_groups(path: &Path) -> Result<Vec<Vec<String>>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut groups: Vec<Vec<String>> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let species: Vec<String> = trimmed
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();
        if species.iter().any(|s| s.is_empty()) {
            return Err(format!(
                "Invalid format in species groups file '{}'",
                path.display()
            )
            .into());
        }
        groups.push(species);
    }
    if groups.is_empty() {
        return Err("Species groups file is empty".into());
    }
    if groups.len() % 2 != 0 {
        return Err(format!(
            "Species groups file '{}' must have an even number of lines",
            path.display()
        )
        .into());
    }
    let mut combos: Vec<Vec<String>> = vec![Vec::new()];
    for group in groups {
        let mut new_combos = Vec::new();
        for combo in &combos {
            for species in &group {
                let mut new_combo = combo.clone();
                new_combo.push(species.clone());
                new_combos.push(new_combo);
            }
        }
        combos = new_combos;
    }
    println!("Generated {} species combinations from file.", combos.len());
    Ok(combos)
}

fn read_genes_list(path: &Path) -> Result<HashSet<String>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            set.insert(trimmed.to_string());
        }
    }
    Ok(set)
}

fn read_fasta(path: &Path) -> Result<(HashMap<String, String>, usize), Box<dyn Error>> {
    // Support both 2-line and multi-line FASTA by concatenating sequence lines
    let content = fs::read_to_string(path)?;
    let mut map: HashMap<String, String> = HashMap::new();
    let mut seq_len: Option<usize> = None;
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            // Flush previous record if any
            if let Some(name) = current_name.take() {
                let cur_len = current_seq.len();
                if let Some(expected) = seq_len {
                    if cur_len != expected {
                        return Err(format!(
                            "Inconsistent sequence length in {}",
                            path.display()
                        )
                        .into());
                    }
                } else {
                    seq_len = Some(cur_len);
                }
                map.insert(name, current_seq.clone());
                current_seq.clear();
            }
            current_name = Some(line[1..].trim().to_string());
        } else {
            current_seq.push_str(line.trim());
        }
    }

    // Flush the last record
    if let Some(name) = current_name.take() {
        let cur_len = current_seq.len();
        if let Some(expected) = seq_len {
            if cur_len != expected {
                return Err(format!(
                    "Inconsistent sequence length in {}",
                    path.display()
                )
                .into());
            }
        } else {
            seq_len = Some(cur_len);
        }
        map.insert(name, current_seq.clone());
    }

    let len = seq_len.unwrap_or(0);
    Ok((map, len))
}

fn write_fasta(path: &Path, species: &[String], seqs: &[String]) -> Result<(), Box<dyn Error>> {
    let mut f = File::create(path)?;
    for (name, seq) in species.iter().zip(seqs.iter()) {
        writeln!(f, ">{name}")?;
        writeln!(f, "{seq}")?;
    }
    Ok(())
}

fn process_alignment(
    map: &HashMap<String, String>,
    seq_len: usize,
    species_list: &[String],
    args: &Args,
) -> Option<(Vec<String>, bool, bool)> {
    let mut sequences: Vec<Vec<char>> = Vec::new();
    let mut species_canceled: Vec<bool> = Vec::new();
    for s in species_list {
        if let Some(seq) = map.get(s) {
            sequences.push(seq.chars().collect());
            species_canceled.push(false);
        } else {
            sequences.push(vec!['-'; seq_len]);
            species_canceled.push(true);
        }
    }

    if args.nix_full_deletions && species_canceled.iter().any(|&c| c) {
        return None;
    }

    let mut fully_canceled_gene = false;
    if species_list.len() % 2 == 0 {
        let mut num_uncanceled = 0;
        for pair in species_canceled.chunks(2) {
            if pair.len() == 2 && !pair[0] && !pair[1] {
                num_uncanceled += 1;
            }
        }
        if num_uncanceled < args.min_pairs {
            fully_canceled_gene = true;
        }
    }

    let outgroup_seq = args
        .outgroup_species
        .as_ref()
        .and_then(|sp| map.get(sp))
        .map(|s| s.chars().collect::<Vec<char>>());

    for idx in 0..seq_len {
        let position_list: Vec<char> = sequences.iter().map(|seq| seq[idx]).collect();
        let cancel_site_due_to_gap = position_list.iter().any(|&c| c == '-');

        if args.cancel_only_partner && cancel_site_due_to_gap {
            let mut pairs_left_after_partner_cancellation = sequences.len() / 2;
            for (pair_idx, pair) in position_list.chunks(2).enumerate() {
                if pair.contains(&'-') {
                    sequences[2 * pair_idx][idx] = '-';
                    sequences[2 * pair_idx + 1][idx] = '-';
                    pairs_left_after_partner_cancellation -= 1;
                }
            }
            if pairs_left_after_partner_cancellation < args.min_pairs {
                for seq in sequences.iter_mut() {
                    seq[idx] = '-';
                }
                continue;
            }
        } else if cancel_site_due_to_gap {
            for seq in sequences.iter_mut() {
                seq[idx] = '-';
            }
            continue;
        }

        if let Some(out_seq) = &outgroup_seq {
            if !cancel_site_due_to_gap {
                let out_res = out_seq[idx];
                if out_res != '-' {
                    let mut mismatch = false;
                    for seq in sequences.iter().skip(1).step_by(2) {
                        if seq[idx] != out_res {
                            mismatch = true;
                            break;
                        }
                    }
                    if mismatch {
                        for seq in sequences.iter_mut() {
                            seq[idx] = '-';
                        }
                        continue;
                    }
                }
            }
        }

        if args.cancel_tri_allelic && species_list.len() == 4 {
            let mut set = HashSet::new();
            for seq in &sequences {
                set.insert(seq[idx]);
            }
            if set.len() == 3 {
                for seq in sequences.iter_mut() {
                    seq[idx] = '-';
                }
            }
        }
    }

    let seq_strings: Vec<String> = sequences
        .iter()
        .map(|chars| chars.iter().collect())
        .collect();
    let gene_has_data = seq_strings
        .iter()
        .any(|s| s.chars().any(|c| c != '-'));

    Some((seq_strings, gene_has_data, fully_canceled_gene))
}

fn run(args: Args) -> Result<(), Box<dyn Error>> {
    if args.response_file.is_none()
        && args.species_groups_file.is_none()
        && args.response_dir.is_none()
    {
        return Err(
            "must provide --response-file, --species-groups-file, or --response-dir".into(),
        );
    }

    let canceled_dir = if let Some(dir) = &args.canceled_alignments_dir {
        dir.clone()
    } else {
        let parent = args
            .alignments_dir
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."));
        let base = if let Some(rf) = &args.response_file {
            rf.file_stem().unwrap().to_string_lossy().to_string()
        } else if let Some(rd) = &args.response_dir {
            rd.file_name().unwrap().to_string_lossy().to_string()
        } else {
            args.species_groups_file
                .as_ref()
                .unwrap()
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .to_string()
        };
        parent.join(format!("{}_gap-canceled_alignments", base))
    };

    fs::create_dir_all(&canceled_dir)?;
    println!("New alignments folder: {}", canceled_dir.display());

    let limited_genes = if let Some(path) = &args.limited_genes_list {
        Some(read_genes_list(path)?)
    } else {
        None
    };

    let species_combos = if let Some(rf) = &args.response_file {
        vec![read_response(rf)?]
    } else if let Some(dir) = &args.response_dir {
        // Deterministic ordering: sort by filename
        let mut entries: Vec<_> = fs::read_dir(dir)?
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
            .collect();
        entries.sort_by_key(|e| e.file_name());

        let mut combos = Vec::new();
        for entry in entries {
            let path = entry.path();
            combos.push(read_response(&path)?);
        }
        combos
    } else {
        parse_species_groups(args.species_groups_file.as_ref().unwrap())?
    };

    for (combo_idx, species_list) in species_combos.iter().enumerate() {
        let combo_dir = if species_combos.len() > 1 {
            let dir = canceled_dir.join(format!("combo_{}-alignments", combo_idx));
            if dir.exists() {
                fs::remove_dir_all(&dir)?;
            }
            fs::create_dir_all(&dir)?;
            dir
        } else {
            canceled_dir.clone()
        };

        println!("Generating alignments for: {}", species_list.join(" "));

        let mut file_count = 0;
        let mut combo_has_valid_gene = false;
        let mut fully_canceled_genes = 0;

        for entry in fs::read_dir(&args.alignments_dir)? {
            let entry = entry?;
            let file_name = entry.file_name();
            let file_name_str = file_name.to_string_lossy().to_string();
            file_count += 1;
            if entry
                .path()
                .extension()
                .and_then(|s| s.to_str())
                != Some("fas")
            {
                continue;
            }
            if let Some(genes) = &limited_genes {
                if !genes.contains(&file_name_str) {
                    continue;
                }
            }

            let (map, seq_len) = match read_fasta(&entry.path()) {
                Ok(res) => res,
                Err(e) => {
                    println!(
                        "WARNING: Skipping file '{}' because it could not be parsed as FASTA. Error: {}",
                        file_name_str, e
                    );
                    continue;
                }
            };
            if map.is_empty() {
                println!("WARNING: Skipping empty alignment file: {}", file_name_str);
                continue;
            }

            if let Some((seqs, gene_has_data, fully_canceled)) =
                process_alignment(&map, seq_len, species_list, &args)
            {
                if gene_has_data {
                    combo_has_valid_gene = true;
                }
                if fully_canceled {
                    fully_canceled_genes += 1;
                }
                let out_path = combo_dir.join(&file_name_str);
                write_fasta(&out_path, species_list, &seqs)?;
            }
        }

        if file_count > 0 && !combo_has_valid_gene {
            return Err(format!(
                "FATAL: For species combo {}, all generated alignment files would consist entirely of gaps ('-').",
                species_list.join(" ")
            )
            .into());
        }
        println!("number of genes fully canceled: {}", fully_canceled_genes);
    }

    println!("finished generating gap-canceled alignments!");
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let start = Instant::now();
    let res = run(args);
    println!("Elapsed: {:.2}s", start.elapsed().as_secs_f64());
    res
}

