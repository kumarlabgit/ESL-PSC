use anyhow::{anyhow, bail, Context, Result};
use chrono::Local;
use clap::Parser;
use csv::Writer;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use rand::Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::time::{Duration, Instant};

fn emit_status(args: std::fmt::Arguments<'_>) {
    let mut stdout = std::io::stdout().lock();
    let _ = stdout.write_fmt(args);
    let _ = stdout.write_all(b"\n");
    let _ = stdout.flush();
}

macro_rules! statusln {
    ($($arg:tt)*) => {{
        emit_status(format_args!($($arg)*));
    }};
}

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(name = "esl-psc")]
#[command(
    about = "Unified in-memory ESL-PSC CLI (preprocess + sparse-group-lasso + multimatrix orchestration)"
)]
struct Args {
    #[arg(long = "input_alignments_dir", alias = "input-alignments-dir")]
    input_alignments_dir: Option<PathBuf>,

    #[arg(long = "alignments_dir", alias = "alignments-dir")]
    alignments_dir: Option<PathBuf>,

    #[arg(
        long = "prediction_alignments_dir",
        alias = "prediction-alignments-dir"
    )]
    prediction_alignments_dir: Option<PathBuf>,

    #[arg(long = "response_matrix_path", alias = "response-matrix-path")]
    response_matrix_path: Option<PathBuf>,

    #[arg(long = "response_file", alias = "response-file")]
    response_file: Option<PathBuf>,

    #[arg(long = "response_dir", alias = "response-dir")]
    response_dir: Option<PathBuf>,

    #[arg(long = "species_groups_file", alias = "species-groups-file")]
    species_groups_file: Option<PathBuf>,

    #[arg(long = "species_pheno_path", alias = "species-pheno-path")]
    species_pheno_path: Option<PathBuf>,

    #[arg(long = "output_dir", alias = "output-dir")]
    output_dir: Option<PathBuf>,

    #[arg(long = "output_file_base_name", alias = "output-file-base-name")]
    output_file_base_name: String,

    #[arg(
        long = "initial_lambda1",
        alias = "initial-lambda1",
        default_value_t = 0.01
    )]
    initial_lambda1: f64,

    #[arg(
        long = "final_lambda1",
        alias = "final-lambda1",
        default_value_t = 0.99
    )]
    final_lambda1: f64,

    #[arg(
        long = "initial_lambda2",
        alias = "initial-lambda2",
        default_value_t = 0.01
    )]
    initial_lambda2: f64,

    #[arg(
        long = "final_lambda2",
        alias = "final-lambda2",
        default_value_t = 0.99
    )]
    final_lambda2: f64,

    #[arg(long = "lambda_step", alias = "lambda-step", default_value_t = 0.05)]
    lambda_step: f64,

    #[arg(long = "use_logspace", alias = "use-logspace", default_value_t = false)]
    use_logspace: bool,

    #[arg(
        long = "num_log_points",
        alias = "num-log-points",
        default_value_t = 20
    )]
    num_log_points: usize,

    #[arg(
        long = "group_penalty_type",
        alias = "group-penalty-type",
        default_value = "median"
    )]
    group_penalty_type: String,

    #[arg(
        long = "initial_gp_value",
        alias = "initial-gp-value",
        default_value_t = 1.0
    )]
    initial_gp_value: f64,

    #[arg(
        long = "final_gp_value",
        alias = "final-gp-value",
        default_value_t = 1.0
    )]
    final_gp_value: f64,

    #[arg(long = "gp_step", alias = "gp-step", default_value_t = 6.0)]
    gp_step: f64,

    #[arg(
        long = "use_continuous_phenotypes",
        alias = "use-continuous-phenotypes",
        default_value_t = false
    )]
    use_continuous_phenotypes: bool,

    #[arg(long = "maxiter", default_value_t = 100)]
    maxiter: usize,

    #[arg(
        long = "disable_ec",
        alias = "disable-ec",
        default_value_t = false,
        help = "Disable epsilon-comparison line-search acceptance in the sparse-group-lasso solver"
    )]
    disable_ec: bool,

    #[arg(long = "num_threads", alias = "num-threads")]
    num_threads: Option<usize>,

    #[arg(
        long = "top_rank_frac",
        alias = "top-rank-frac",
        default_value_t = 0.01
    )]
    top_rank_frac: f64,

    #[arg(
        long = "no_pred_output",
        alias = "no-pred-output",
        default_value_t = false
    )]
    no_pred_output: bool,

    #[arg(
        long = "no_genes_output",
        alias = "no-genes-output",
        default_value_t = false
    )]
    no_genes_output: bool,

    #[arg(
        long = "show_selected_sites",
        alias = "show-selected-sites",
        default_value_t = false
    )]
    show_selected_sites: bool,

    #[arg(
        long = "use_uncanceled_alignments",
        alias = "use-uncanceled-alignments",
        default_value_t = false
    )]
    use_uncanceled_alignments: bool,

    #[arg(
        long = "cancel_tri_allelic",
        alias = "cancel-tri-allelic",
        default_value_t = false
    )]
    cancel_tri_allelic: bool,

    #[arg(
        long = "nix_full_deletions",
        alias = "nix-full-deletions",
        default_value_t = false
    )]
    nix_full_deletions: bool,

    #[arg(long = "outgroup_species", alias = "outgroup-species")]
    outgroup_species: Option<String>,

    #[arg(
        long = "cancel_only_partner",
        alias = "cancel-only-partner",
        default_value_t = false
    )]
    cancel_only_partner: bool,

    #[arg(long = "min_pairs", alias = "min-pairs", default_value_t = 2)]
    min_pairs: usize,

    #[arg(long = "limited_genes_list", alias = "limited-genes-list")]
    limited_genes_list: Option<PathBuf>,

    // Compatibility placeholders for full ESL-PSC CLI surface.
    #[arg(long = "esl_inputs_outputs_dir", alias = "esl-inputs-outputs-dir")]
    esl_inputs_outputs_dir: Option<PathBuf>,

    #[arg(long = "esl_main_dir", alias = "esl-main-dir")]
    esl_main_dir: Option<PathBuf>,

    #[arg(long = "path_file_path", alias = "path-file-path")]
    path_file_path: Option<PathBuf>,

    #[arg(long = "preprocessed_dir_name", alias = "preprocessed-dir-name")]
    preprocessed_dir_name: Option<String>,

    #[arg(long = "pheno_names", alias = "pheno-names", num_args = 2)]
    pheno_names: Option<Vec<String>>,

    #[arg(long = "min_genes", alias = "min-genes", default_value_t = 0)]
    min_genes: usize,

    #[arg(
        long = "use_existing_preprocess",
        alias = "use-existing-preprocess",
        default_value_t = false
    )]
    use_existing_preprocess: bool,

    #[arg(
        long = "use_default_gp",
        alias = "use-default-gp",
        default_value_t = false
    )]
    use_default_gp: bool,

    #[arg(
        long = "keep_raw_output",
        alias = "keep-raw-output",
        default_value_t = false
    )]
    keep_raw_output: bool,

    #[arg(
        long = "make_sps_plot",
        alias = "make-sps-plot",
        default_value_t = false
    )]
    make_sps_plot: bool,

    #[arg(
        long = "make_sps_kde_plot",
        alias = "make-sps-kde-plot",
        default_value_t = false
    )]
    make_sps_kde_plot: bool,

    #[arg(
        long = "make_continuous_plot",
        alias = "make-continuous-plot",
        default_value_t = false
    )]
    make_continuous_plot: bool,

    #[arg(long = "canceled_alignments_dir", alias = "canceled-alignments-dir")]
    canceled_alignments_dir: Option<PathBuf>,

    #[arg(
        long = "use_existing_alignments",
        alias = "use-existing-alignments",
        default_value_t = false
    )]
    use_existing_alignments: bool,

    #[arg(
        long = "delete_preprocess",
        alias = "delete-preprocess",
        default_value_t = false
    )]
    delete_preprocess: bool,

    #[arg(
        long = "preserve_canceled_alignments",
        alias = "preserve-canceled-alignments",
        default_value_t = false
    )]
    preserve_canceled_alignments: bool,

    #[arg(
        long = "make_null_models",
        alias = "make-null-models",
        default_value_t = false
    )]
    make_null_models: bool,

    #[arg(
        long = "make_pair_randomized_null_models",
        alias = "make-pair-randomized-null-models",
        default_value_t = false
    )]
    make_pair_randomized_null_models: bool,

    #[arg(
        long = "num_randomized_alignments",
        alias = "num-randomized-alignments",
        default_value_t = 10
    )]
    num_randomized_alignments: usize,

    #[arg(
        long = "auto_convert_to_2line",
        alias = "auto-convert-to-2line",
        default_value_t = false
    )]
    auto_convert_to_2line: bool,

    #[arg(
        long = "no_checkpoint",
        alias = "no-checkpoint",
        default_value_t = false
    )]
    no_checkpoint: bool,

    #[arg(
        long = "force_from_beginning",
        alias = "force-from-beginning",
        default_value_t = false
    )]
    force_from_beginning: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FeatureMeta {
    label: String,
    gene_idx: usize,
    position: usize,
    aa: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GeneMeta {
    name: String,
    feature_start: usize,
    feature_end: usize,
    var_site_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PreprocessedData {
    species: Vec<String>,
    y_model: Vec<f64>,
    x: Vec<Vec<f64>>, // shape: n_samples x n_features
    features: Vec<FeatureMeta>,
    genes: Vec<GeneMeta>,
}

#[derive(Debug)]
struct GeneFeatureColumn {
    position: usize,
    aa: u8,
    values: Vec<f64>,
}

#[derive(Debug)]
struct GenePreprocessResult {
    gene_name: String,
    var_site_count: usize,
    columns: Vec<GeneFeatureColumn>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ModelResult {
    lambda1: f64,
    lambda2: f64,
    penalty_term: f64,
    beta: Vec<f64>,
    intercept: f64,
    rmse: f64,
    num_genes: usize,
    gene_gss: Vec<f64>,
    selected_sites: Vec<SelectedSite>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SelectedSite {
    gene_idx: usize,
    position: usize,
    pss: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GeneAlignment {
    name: String,
    seq_len: usize,
    seqs: HashMap<String, Vec<u8>>, // species -> uppercase sequence
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ComboJob {
    index: usize,
    combo_label: String,
    combo_tag: String,
    species: Vec<String>,
    y_raw: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PhenotypeInfo {
    values: HashMap<String, f64>,
    is_binary: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PredictionDesign {
    species: Vec<String>,
    // For each feature index, row indices (within `species`) with a match.
    feature_hit_rows: Vec<Vec<usize>>,
    true_values: Vec<Option<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PredictionRowOut {
    species_combo: String,
    lambda1: f64,
    lambda2: f64,
    penalty_term: f64,
    num_genes: usize,
    input_rmse: f64,
    species: String,
    sps: f64,
    true_phenotype: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GeneAggregate {
    name: String,
    single_highest_gss: f64,
    single_best_rank: Option<usize>,
    highest_ever_gss: f64,
    best_ever_rank: Option<usize>,
    num_combos_ranked: usize,
    num_combos_ranked_top: usize,
    selected_sites: HashMap<usize, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CheckpointState {
    next_combo_index: usize,
    total_model_runs: usize,
    gene_aggregates: Vec<GeneAggregate>,
    prediction_rows: Vec<PredictionRowOut>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CheckpointRunAudit {
    combo: usize,
    lambda1: f64,
    lambda2: f64,
    penalty_term: f64,
    input_rmse: f64,
}

fn main() -> Result<()> {
    let start = Instant::now();
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(subcommand) = raw_args.get(1).map(|s| s.as_str()) {
        match subcommand {
            "pairs" => {
                return run_python_module_subcommand(
                    "esl_psc_cli.auto_pairs_cli",
                    &raw_args[2..],
                    &[],
                )
            }
            "plot" => {
                return run_python_module_subcommand("esl_psc_cli.plot_cli", &raw_args[2..], &[])
            }
            "site-counter" => {
                let current_exe = std::env::current_exe()
                    .context("unable to determine path to current esl-psc executable")?;
                let envs = vec![
                    ("SITE_COUNTER_RS", current_exe.to_string_lossy().to_string()),
                    ("SITE_COUNTER_RS_DISABLE", "0".to_string()),
                ];
                return run_python_module_subcommand(
                    "esl_psc_cli.fast_scan_cli",
                    &raw_args[2..],
                    &envs,
                );
            }
            "site-counter-backend" => {
                return site_counter_rs::run_backend_stdio().map_err(|e| anyhow!(e));
            }
            "run" => {
                let args = parse_args_with_config(&raw_args[2..]);
                return run_unified_pipeline(args, start);
            }
            _ => {}
        }
    }

    let args = parse_args_with_config(raw_args.get(1..).unwrap_or(&[]));
    run_unified_pipeline(args, start)
}

fn run_unified_pipeline(args: Args, start: Instant) -> Result<()> {
    let base_alignments = resolve_alignments_dir(&args)?;
    let base_prediction = resolve_prediction_alignments_dir(&args, &base_alignments);
    let (resolved_alignments, resolved_prediction) = if args.auto_convert_to_2line {
        let converted_alignments = ensure_two_line_dir(&base_alignments, false)?;
        let converted_prediction = if base_prediction == base_alignments {
            converted_alignments.clone()
        } else {
            ensure_two_line_dir(&base_prediction, false)?
        };
        (converted_alignments, converted_prediction)
    } else {
        (base_alignments, base_prediction)
    };
    validate_args(&args, &resolved_alignments)?;
    report_compat_warnings(&args);
    let lambda_grid = build_lambda_grid(&args)?;
    let model_workers = configure_model_workers(&args, lambda_grid.len())?;

    let output_dir = resolve_output_dir(&args)?;
    fs::create_dir_all(&output_dir)
        .with_context(|| format!("failed to create output dir {}", output_dir.display()))?;
    let preprocess_root = args
        .esl_inputs_outputs_dir
        .clone()
        .unwrap_or_else(|| output_dir.join("preprocessed_data_and_models"));
    let use_preprocess_dirs = args.use_existing_preprocess
        || args.esl_inputs_outputs_dir.is_some()
        || args.preprocessed_dir_name.is_some()
        || args.delete_preprocess
        || args.keep_raw_output;
    if use_preprocess_dirs {
        fs::create_dir_all(&preprocess_root).with_context(|| {
            format!(
                "failed to create preprocess dir root {}",
                preprocess_root.display()
            )
        })?;
    }
    let model_dir = preprocess_root.clone();
    write_run_config_txt(&output_dir, &args)?;

    let limited_gene_set = if let Some(path) = &args.limited_genes_list {
        Some(read_limited_gene_set(path)?)
    } else {
        None
    };

    statusln!(
        "Loading input alignments from {}",
        resolved_alignments.display()
    );
    let path_file_order = if let Some(path_file) = &args.path_file_path {
        Some(read_path_file_order(path_file, &resolved_alignments)?)
    } else {
        None
    };

    let train_alignments = load_alignments(
        &resolved_alignments,
        limited_gene_set.as_ref(),
        path_file_order.as_deref(),
        Some("input"),
    )?;
    if train_alignments.is_empty() {
        if args.use_existing_alignments {
            statusln!(
                "No FASTA files found in base alignments directory; using per-combo existing alignments."
            );
        } else {
            bail!(
                "no usable alignments found in {}",
                resolved_alignments.display()
            );
        }
    } else {
        statusln!("Loaded {} gene alignments", train_alignments.len());
    }

    let prediction_alignments_owned: Option<Vec<GeneAlignment>> =
        if args.no_pred_output || resolved_prediction == resolved_alignments {
            None
        } else {
            statusln!(
                "Loading prediction alignments from {}",
                resolved_prediction.display()
            );
            Some(load_alignments(
                &resolved_prediction,
                limited_gene_set.as_ref(),
                None,
                Some("prediction"),
            )?)
        };
    let prediction_alignments: &[GeneAlignment] = if args.no_pred_output {
        &[]
    } else if resolved_prediction == resolved_alignments {
        // Reuse loaded input alignments directly when prediction and input
        // directories are the same to avoid duplicating all sequence data.
        train_alignments.as_slice()
    } else {
        prediction_alignments_owned
            .as_deref()
            .ok_or_else(|| anyhow!("failed to load prediction alignments"))?
    };

    let phenotype_info = if let Some(path) = &args.species_pheno_path {
        Some(read_species_phenotypes(path)?)
    } else {
        None
    };

    let mut combos = build_combo_jobs(&args, phenotype_info.as_ref())?;
    if args.make_null_models {
        combos = make_null_combo_jobs(&combos)?;
    }
    if combos.is_empty() {
        bail!("no combo jobs were created");
    }

    if !train_alignments.is_empty() {
        validate_combo_species_against_alignments(&combos, &train_alignments)?;
    }

    let is_multimatrix_mode = args.response_dir.is_some() || args.species_groups_file.is_some();
    let apply_gap_cancel = is_multimatrix_mode && !args.use_uncanceled_alignments;
    let any_continuous_response = combos
        .iter()
        .flat_map(|c| c.y_raw.iter().copied())
        .any(|v| !is_binary_pheno_value(v));
    let auto_continuous = args.response_dir.is_some() && any_continuous_response;
    let use_continuous = args.use_continuous_phenotypes || auto_continuous;

    if auto_continuous && !args.use_continuous_phenotypes {
        statusln!("Detected non-binary response values; using continuous solver mode automatically");
    }
    let preserve_canceled_root = if args.preserve_canceled_alignments
        && is_multimatrix_mode
        && !args.use_existing_alignments
        && !args.make_pair_randomized_null_models
    {
        let dir = derive_canceled_alignments_dir(&args, &output_dir);
        fs::create_dir_all(&dir).with_context(|| {
            format!("failed to create canceled alignments dir {}", dir.display())
        })?;
        Some(dir)
    } else {
        None
    };

    statusln!("Lambda grid has {} pairs", lambda_grid.len());
    statusln!("Lambda model workers: {}", model_workers);

    let mut gene_aggregates: Vec<GeneAggregate> = train_alignments
        .iter()
        .map(|g| GeneAggregate {
            name: g.name.clone(),
            single_highest_gss: 0.0,
            single_best_rank: None,
            highest_ever_gss: 0.0,
            best_ever_rank: None,
            num_combos_ranked: 0,
            num_combos_ranked_top: 0,
            selected_sites: HashMap::new(),
        })
        .collect();

    let mut all_prediction_rows: Vec<PredictionRowOut> = Vec::new();
    let mut top_rank_threshold = (gene_aggregates.len() as f64 * args.top_rank_frac).max(1.0);

    let mut total_model_runs = 0usize;
    let mut start_combo_index = 0usize;
    let checkpoint_enabled = !args.no_checkpoint && combos.len() > 1;
    let checkpoint_min_interval = Duration::from_secs(60);
    let mut checkpoint_last_save = Instant::now();
    let mut pending_checkpoint_audits: Vec<CheckpointRunAudit> = Vec::new();
    if checkpoint_enabled {
        if let Some(restored) = restore_checkpoint_if_available(
            &args,
            &output_dir,
            &mut gene_aggregates,
            &mut all_prediction_rows,
            &mut total_model_runs,
        )? {
            start_combo_index = restored;
            if start_combo_index >= combos.len() {
                statusln!("All combos completed according to checkpoint; skipping integration.");
            } else if start_combo_index > 0 {
                statusln!(
                    "Resuming from checkpoint at combo {} of {}",
                    start_combo_index + 1,
                    combos.len()
                );
            }
        }
    }

    let mut combo_alignment_cache: HashMap<PathBuf, Vec<GeneAlignment>> = HashMap::new();
    for combo in combos.iter().skip(start_combo_index) {
        statusln!(
            "\n--- Processing combo {} of {} ({}) ---",
            combo.index + 1,
            combos.len(),
            combo.combo_label
        );

        let combo_alignments = resolve_combo_alignments(
            &args,
            combo,
            &train_alignments,
            limited_gene_set.as_ref(),
            &mut combo_alignment_cache,
        )?;
        validate_combo_species_against_alignments(std::slice::from_ref(combo), combo_alignments)?;
        if gene_aggregates.is_empty() {
            gene_aggregates = combo_alignments
                .iter()
                .map(|g| GeneAggregate {
                    name: g.name.clone(),
                    single_highest_gss: 0.0,
                    single_best_rank: None,
                    highest_ever_gss: 0.0,
                    best_ever_rank: None,
                    num_combos_ranked: 0,
                    num_combos_ranked_top: 0,
                    selected_sites: HashMap::new(),
                })
                .collect();
            top_rank_threshold = (gene_aggregates.len() as f64 * args.top_rank_frac).max(1.0);
        }
        if combo_alignments.len() != gene_aggregates.len() {
            bail!(
                "gene count mismatch across combos (expected {}, got {})",
                gene_aggregates.len(),
                combo_alignments.len()
            );
        }
        let mut combo_highest_gss = vec![0.0_f64; gene_aggregates.len()];
        let mut combo_best_rank: Vec<Option<usize>> = vec![None; gene_aggregates.len()];
        let combo_preserve_dir = preserve_canceled_root.as_ref().map(|root| {
            if combos.len() > 1 {
                root.join(format!("{}-alignments", combo.combo_label))
            } else {
                root.clone()
            }
        });
        if let Some(dir) = &combo_preserve_dir {
            if dir.exists() {
                fs::remove_dir_all(dir).with_context(|| {
                    format!("failed to clear canceled alignments dir {}", dir.display())
                })?;
            }
            fs::create_dir_all(dir).with_context(|| {
                format!("failed to create canceled alignments dir {}", dir.display())
            })?;
        }

        let random_repeats = if args.make_pair_randomized_null_models {
            args.num_randomized_alignments.max(1)
        } else {
            1
        };
        let mut combo_run_audits: Vec<CheckpointRunAudit> = Vec::new();

        for random_rep in 0..random_repeats {
            let preprocess_name = preprocess_dir_name_for_combo(
                &args,
                combo,
                &resolved_alignments,
                is_multimatrix_mode,
            );
            let preprocess_dir = preprocess_root.join(&preprocess_name);
            let prep = if args.use_existing_preprocess {
                load_preprocessed_data(&preprocess_dir, combo, combo_alignments, use_continuous)?
            } else {
                let p = preprocess_combo_alignments(
                    combo_alignments,
                    combo,
                    use_continuous,
                    apply_gap_cancel && !args.use_existing_alignments,
                    args.make_pair_randomized_null_models,
                    combo_preserve_dir.as_deref(),
                    &args,
                )?;
                if use_preprocess_dirs {
                    write_preprocessed_data(&preprocess_dir, &p)?;
                }
                p
            };

            if prep.genes.len() != gene_aggregates.len() {
                bail!(
                    "gene count mismatch between preprocessed data and combo alignments (expected {}, got {})",
                    gene_aggregates.len(),
                    prep.genes.len()
                );
            }

            let rep_label = if args.make_pair_randomized_null_models {
                format!("{}_{}", combo.combo_label, random_rep)
            } else {
                combo.combo_label.clone()
            };
            maybe_dump_preprocess_features(&output_dir, &rep_label, &prep)?;

            let penalty_terms = if args.group_penalty_type.eq_ignore_ascii_case("median") {
                if apply_gap_cancel && !args.use_existing_alignments {
                    build_penalty_terms(&args, &prep.genes)?
                } else {
                    vec![median_var_sites_from_alignments(combo_alignments)?]
                }
            } else {
                build_penalty_terms(&args, &prep.genes)?
            };
            statusln!("Penalty schedule for {}: {:?}", rep_label, penalty_terms);

            let pred_design = if args.no_pred_output {
                None
            } else {
                Some(build_prediction_design(
                    prediction_alignments,
                    &prep.features,
                    &prep.genes,
                    &combo.species,
                    phenotype_info.as_ref(),
                )?)
            };
            let input_feature_hit_rows = if args.no_pred_output {
                None
            } else {
                Some(build_feature_hit_rows_from_dense(&prep.x))
            };

            let lipschitz = if prep.features.is_empty() {
                1.0
            } else {
                estimate_lipschitz(&prep.x, use_continuous)
            };
            let lambda_rows = lambda_row_ranges(&lambda_grid);

            for penalty in &penalty_terms {
                let effective_group_penalty_kind = if args.use_default_gp {
                    "std"
                } else {
                    args.group_penalty_type.as_str()
                };
                let group_weights =
                    compute_group_weights(effective_group_penalty_kind, *penalty, &prep.genes);
                let combo_tag = if args.make_pair_randomized_null_models {
                    format!("{}_{}", combo.combo_tag, random_rep)
                } else {
                    combo.combo_tag.clone()
                };
                statusln!("Building models...");
                let grid_run_total = lambda_grid.len();
                let progress_counter = Arc::new(AtomicUsize::new(0));
                let progress_done = Arc::new(AtomicBool::new(false));
                let progress_counter_for_thread = Arc::clone(&progress_counter);
                let progress_done_for_thread = Arc::clone(&progress_done);
                let progress_handle = std::thread::spawn(move || {
                    let mut last_reported = 0usize;
                    loop {
                        let current = progress_counter_for_thread
                            .load(AtomicOrdering::Relaxed)
                            .min(grid_run_total);
                        if current != last_reported {
                            statusln!(
                                "run {} of {} in current grid;  time: {}",
                                current,
                                grid_run_total,
                                Local::now().format("%H:%M:%S")
                            );
                            last_reported = current;
                        }
                        if progress_done_for_thread.load(AtomicOrdering::Relaxed) {
                            if last_reported < grid_run_total {
                                statusln!(
                                    "run {} of {} in current grid;  time: {}",
                                    grid_run_total,
                                    grid_run_total,
                                    Local::now().format("%H:%M:%S")
                                );
                            }
                            break;
                        }
                        std::thread::sleep(Duration::from_millis(100));
                    }
                });

                // Solve lambda grids by lambda1 rows while parallelizing rows.
                let row_results: Vec<Vec<ModelResult>> =
                    if lambda_rows.len() <= 1 || model_workers <= 1 {
                        lambda_rows
                            .iter()
                            .map(|row| {
                                solve_lambda_row(
                                    &prep,
                                    &group_weights,
                                    &lambda_grid[row.clone()],
                                    *penalty,
                                    use_continuous,
                                    args.maxiter,
                                    args.disable_ec,
                                    lipschitz,
                                    Some(Arc::clone(&progress_counter)),
                                )
                            })
                            .collect::<Result<Vec<_>>>()?
                    } else {
                        lambda_rows
                            .par_iter()
                            .map(|row| {
                                solve_lambda_row(
                                    &prep,
                                    &group_weights,
                                    &lambda_grid[row.clone()],
                                    *penalty,
                                    use_continuous,
                                    args.maxiter,
                                    args.disable_ec,
                                    lipschitz,
                                    Some(Arc::clone(&progress_counter)),
                                )
                            })
                            .collect::<Result<Vec<_>>>()?
                    };
                progress_done.store(true, AtomicOrdering::Relaxed);
                if let Err(e) = progress_handle.join() {
                    bail!("model-progress thread failed: {:?}", e);
                }

                for results in row_results {
                    for result in &results {
                        if args.keep_raw_output {
                            write_model_file(
                                &model_dir,
                                &args.output_file_base_name,
                                &prep.features,
                                result,
                                if is_multimatrix_mode {
                                    Some(rep_label.as_str())
                                } else {
                                    None
                                },
                            )?;
                        }

                        update_gene_stats_for_run(
                            result,
                            &mut combo_highest_gss,
                            &mut combo_best_rank,
                            &mut gene_aggregates,
                        );

                        if let Some(pred) = &pred_design {
                            append_prediction_rows(
                                &mut all_prediction_rows,
                                &prep,
                                input_feature_hit_rows
                                    .as_deref()
                                    .ok_or_else(|| anyhow!("missing input feature hit rows"))?,
                                pred,
                                result,
                                &combo_tag,
                                phenotype_info.as_ref(),
                            );
                        }

                        total_model_runs += 1;
                        combo_run_audits.push(CheckpointRunAudit {
                            combo: combo.index,
                            lambda1: result.lambda1,
                            lambda2: result.lambda2,
                            penalty_term: result.penalty_term,
                            input_rmse: result.rmse,
                        });
                    }
                }
            }
        }

        if is_multimatrix_mode {
            finalize_combo_multimatrix_stats(
                &combo_highest_gss,
                &combo_best_rank,
                top_rank_threshold,
                &mut gene_aggregates,
            );
        }
        if use_preprocess_dirs && args.delete_preprocess && preprocess_root.exists() {
            let preprocess_name = preprocess_dir_name_for_combo(
                &args,
                combo,
                &resolved_alignments,
                is_multimatrix_mode,
            );
            let preprocess_dir = preprocess_root.join(preprocess_name);
            if preprocess_dir.is_dir() {
                fs::remove_dir_all(&preprocess_dir).with_context(|| {
                    format!(
                        "failed deleting preprocess dir {}",
                        preprocess_dir.display()
                    )
                })?;
            }
        }

        if checkpoint_enabled {
            pending_checkpoint_audits.extend(combo_run_audits.into_iter());
            let is_last_combo = combo.index + 1 >= combos.len();
            let save_due = checkpoint_last_save.elapsed() >= checkpoint_min_interval;
            if is_last_combo || save_due {
                save_checkpoint_state(
                    &args,
                    &output_dir,
                    combo.index + 1,
                    total_model_runs,
                    &gene_aggregates,
                    &all_prediction_rows,
                    &pending_checkpoint_audits,
                )?;
                pending_checkpoint_audits.clear();
                checkpoint_last_save = Instant::now();
            }
        }
    }

    statusln!("\nFinished {} model runs", total_model_runs);

    if !args.no_pred_output {
        let pred_path = write_predictions_csv(
            &output_dir,
            &args.output_file_base_name,
            &all_prediction_rows,
            phenotype_info.is_some(),
        )?;
        maybe_generate_plots(&args, &pred_path, use_continuous)?;
    }

    if !args.no_genes_output {
        if is_multimatrix_mode {
            write_gene_ranks_csv_multimatrix(
                &output_dir,
                &args.output_file_base_name,
                &gene_aggregates,
                args.show_selected_sites,
            )?;
        } else {
            write_gene_ranks_csv_single(
                &output_dir,
                &args.output_file_base_name,
                &gene_aggregates,
                args.show_selected_sites,
            )?;
        }

        if args.show_selected_sites {
            write_selected_sites_csv(&output_dir, &args.output_file_base_name, &gene_aggregates)?;
        }
    }

    statusln!(
        "Unified run complete in {:.3}s",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn parse_args_with_config(cli_tail: &[String]) -> Args {
    let cfg = Path::new("esl_psc_config.txt");
    if !cfg.is_file() {
        let mut argv = vec!["esl-psc".to_string()];
        argv.extend(cli_tail.iter().cloned());
        return Args::parse_from(argv);
    }
    statusln!("getting args from esl_psc_config.txt...");
    let cfg_text = fs::read_to_string(cfg).unwrap_or_default();
    let mut merged: Vec<String> = Vec::with_capacity(1 + cli_tail.len() + cfg_text.len() / 4);
    merged.push("esl-psc".to_string());
    merged.extend(cfg_text.split_whitespace().map(|s| s.to_string()));
    merged.extend(cli_tail.iter().cloned());
    Args::parse_from(merged)
}

fn run_python_module_subcommand(
    module_name: &str,
    module_args: &[String],
    extra_env: &[(&str, String)],
) -> Result<()> {
    let mut last_error: Option<String> = None;
    for python_cmd in python_command_candidates() {
        let mut cmd = Command::new(&python_cmd);
        cmd.arg("-m").arg(module_name);
        cmd.args(module_args);
        configure_python_env_for_toolkit(&mut cmd);
        for (key, value) in extra_env {
            cmd.env(key, value);
        }
        match cmd.status() {
            Ok(status) if status.success() => return Ok(()),
            Ok(status) => {
                last_error = Some(format!(
                    "python module '{}' via '{}' exited with status {}",
                    module_name, python_cmd, status
                ));
            }
            Err(e) => {
                last_error = Some(format!(
                    "failed to launch python '{}' for module '{}': {}",
                    python_cmd, module_name, e
                ));
            }
        }
    }

    bail!(
        "{}",
        last_error.unwrap_or_else(|| {
            format!(
                "could not execute python module '{}' (no working python interpreter found)",
                module_name
            )
        })
    );
}

fn configure_python_env_for_toolkit(cmd: &mut Command) {
    let candidate_paths = toolkit_pythonpath_candidates();
    if candidate_paths.is_empty() {
        return;
    }
    let sep = if cfg!(windows) { ";" } else { ":" };
    let existing = std::env::var("PYTHONPATH").unwrap_or_default();
    let mut pieces: Vec<String> = Vec::new();
    for path in candidate_paths {
        pieces.push(path.to_string_lossy().to_string());
    }
    if !existing.trim().is_empty() {
        pieces.push(existing);
    }
    cmd.env("PYTHONPATH", pieces.join(sep));
}

fn toolkit_pythonpath_candidates() -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    if let Ok(path) = std::env::var("ESL_PSC_PYTHONPATH") {
        let trimmed = path.trim();
        if !trimmed.is_empty() {
            let from_env = PathBuf::from(trimmed);
            if from_env.is_dir() {
                candidates.push(from_env);
            }
        }
    }
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let direct = exe_dir.join("python");
            if direct.is_dir() {
                candidates.push(direct);
            }
            if let Some(parent) = exe_dir.parent() {
                let sibling = parent.join("python");
                if sibling.is_dir() {
                    candidates.push(sibling);
                }
                let lib_layout = parent.join("lib").join("esl-psc").join("python");
                if lib_layout.is_dir() {
                    candidates.push(lib_layout);
                }
            }
        }
    }
    dedupe_pathbufs(candidates)
}

fn resolve_alignments_dir(args: &Args) -> Result<PathBuf> {
    if let Some(p) = &args.alignments_dir {
        return Ok(p.clone());
    }
    if let Some(p) = &args.input_alignments_dir {
        return Ok(p.clone());
    }
    if let Some(p) = &args.prediction_alignments_dir {
        return Ok(p.clone());
    }
    if args.use_existing_alignments && args.no_pred_output {
        if let Some(p) = &args.canceled_alignments_dir {
            return Ok(p.clone());
        }
        bail!(
            "--use_existing_alignments with --no_pred_output requires --canceled_alignments_dir when no alignments dir is provided"
        );
    }
    bail!(
        "one of --alignments_dir, --input_alignments_dir, or --prediction_alignments_dir must be provided"
    )
}

fn resolve_prediction_alignments_dir(args: &Args, alignments_dir: &Path) -> PathBuf {
    if let Some(p) = &args.prediction_alignments_dir {
        return p.clone();
    }
    alignments_dir.to_path_buf()
}

fn resolve_output_dir(args: &Args) -> Result<PathBuf> {
    if let Some(p) = &args.output_dir {
        return Ok(p.clone());
    }
    let esl_main_dir = if let Some(p) = &args.esl_main_dir {
        p.clone()
    } else {
        std::env::current_dir().context("unable to get current working directory")?
    };
    let parent = esl_main_dir
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or(esl_main_dir.clone());
    let timestamp = Local::now().format("%Y%m%d_%H%M%S").to_string();
    Ok(parent.join(format!("{}_{}", args.output_file_base_name, timestamp)))
}

fn write_run_config_txt(output_dir: &Path, args: &Args) -> Result<()> {
    let path = output_dir.join(format!("{}_run_config.txt", args.output_file_base_name));
    let value = serde_json::to_value(args)?;
    let obj = value
        .as_object()
        .ok_or_else(|| anyhow!("failed to serialize args object"))?;
    let mut keys: Vec<&str> = obj.keys().map(String::as_str).collect();
    keys.sort_unstable();

    let mut out = String::new();
    for key in keys {
        if key.starts_with('_') {
            continue;
        }
        let flag = format!("--{key}");
        let Some(v) = obj.get(key) else {
            continue;
        };
        match v {
            serde_json::Value::Bool(true) => {
                out.push_str(&flag);
                out.push('\n');
            }
            serde_json::Value::Bool(false) | serde_json::Value::Null => {}
            serde_json::Value::String(s) => {
                out.push_str(&flag);
                out.push(' ');
                out.push_str(s);
                out.push('\n');
            }
            serde_json::Value::Number(n) => {
                out.push_str(&flag);
                out.push(' ');
                out.push_str(&n.to_string());
                out.push('\n');
            }
            serde_json::Value::Array(arr) => {
                if arr.is_empty() {
                    continue;
                }
                out.push_str(&flag);
                for item in arr {
                    out.push(' ');
                    match item {
                        serde_json::Value::String(s) => out.push_str(s),
                        _ => out.push_str(&item.to_string()),
                    }
                }
                out.push('\n');
            }
            _ => {}
        }
    }
    fs::write(&path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn derive_canceled_alignments_dir(args: &Args, output_dir: &Path) -> PathBuf {
    if let Some(p) = &args.canceled_alignments_dir {
        return p.clone();
    }
    let base = if let Some(path) = &args.species_groups_file {
        path.file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("species_groups")
            .replace(".txt", "")
    } else {
        args.output_file_base_name.replace(".txt", "")
    };
    output_dir.join(format!("{}_gap-canceled_alignments", base))
}

fn configure_model_workers(args: &Args, lambda_pairs: usize) -> Result<usize> {
    let n = if let Some(n) = args.num_threads {
        if n == 0 {
            bail!("--num_threads must be >= 1");
        }
        n
    } else {
        auto_model_workers(lambda_pairs)
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .map_err(|e| anyhow!("failed to configure rayon thread pool: {e}"))?;
    Ok(rayon::current_num_threads())
}

fn auto_model_workers(lambda_pairs: usize) -> usize {
    let avail = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    // Keep some headroom to reduce contention with other processes and memory bandwidth pressure.
    let reserve = if avail >= 16 {
        4
    } else if avail >= 8 {
        2
    } else if avail >= 2 {
        1
    } else {
        0
    };
    let mut workers = avail.saturating_sub(reserve).max(1);
    // Never launch more workers than independent lambda models in a combo.
    if lambda_pairs > 0 {
        workers = workers.min(lambda_pairs);
    }
    workers
}

fn validate_args(args: &Args, alignments_dir: &Path) -> Result<()> {
    if !alignments_dir.is_dir() {
        bail!(
            "alignments directory does not exist: {}",
            alignments_dir.display()
        );
    }

    let mut sources = 0usize;
    if args.response_matrix_path.is_some() {
        sources += 1;
    }
    if args.response_file.is_some() {
        sources += 1;
    }
    if args.response_dir.is_some() {
        sources += 1;
    }
    if args.species_groups_file.is_some() {
        sources += 1;
    }
    if sources != 1 {
        bail!(
            "provide exactly one of --response_matrix_path/--response_file, --response_dir, or --species_groups_file"
        );
    }

    if args.initial_lambda1 <= 0.0 || args.final_lambda1 <= 0.0 {
        bail!("lambda1 bounds must be > 0");
    }
    if args.initial_lambda2 <= 0.0 || args.final_lambda2 <= 0.0 {
        bail!("lambda2 bounds must be > 0");
    }
    if args.use_logspace {
        if args.num_log_points < 2 {
            bail!("--num_log_points must be >= 2 for logspace mode");
        }
    } else if args.lambda_step <= 0.0 {
        bail!("--lambda_step must be > 0 for linear grid mode");
    }
    if args.gp_step <= 0.0 {
        bail!("--gp_step must be > 0");
    }
    if let Some(n) = args.num_threads {
        if n == 0 {
            bail!("--num_threads must be >= 1");
        }
    }
    if args.maxiter == 0 {
        bail!("--maxiter must be > 0");
    }
    if args.show_selected_sites && args.no_genes_output {
        bail!("--show_selected_sites cannot be used with --no_genes_output");
    }
    let plot_flags = args.make_sps_plot || args.make_sps_kde_plot || args.make_continuous_plot;
    if plot_flags && args.no_pred_output {
        bail!("plot options require predictions output (disable --no_pred_output)");
    }
    if plot_flags && args.species_pheno_path.is_none() {
        bail!("plot options require --species_pheno_path");
    }
    if args.make_sps_plot && args.make_sps_kde_plot {
        bail!("--make_sps_plot and --make_sps_kde_plot are mutually exclusive");
    }
    if args.num_randomized_alignments == 0 {
        bail!("--num_randomized_alignments must be >= 1");
    }

    let kind = args.group_penalty_type.to_lowercase();
    if !["median", "linear", "sqrt", "std"].contains(&kind.as_str()) {
        bail!("--group_penalty_type must be one of: median, linear, sqrt, std");
    }

    if args.min_pairs == 0 {
        bail!("--min_pairs must be >= 1");
    }
    if let Some(path_file) = &args.path_file_path {
        if !path_file.is_file() {
            bail!("--path_file_path not found: {}", path_file.display());
        }
    }

    Ok(())
}

fn report_compat_warnings(args: &Args) {
    if args.disable_ec {
        statusln!(
            "Compatibility mode: --disable_ec uses strict line-search acceptance. On Linux this reproduces the original ESL-PSC paper-era solver behavior."
        );
    }
}

fn make_null_combo_jobs(combos: &[ComboJob]) -> Result<Vec<ComboJob>> {
    let mut out: Vec<ComboJob> = Vec::new();
    for (combo_num, combo) in combos.iter().enumerate() {
        if combo.species.len() % 4 != 0 {
            bail!("--make_null_models requires an even number of pairs per combo");
        }
        let num_pairs = combo.species.len() / 2;
        let pair_indices: Vec<usize> = (0..num_pairs).collect();
        let k = num_pairs / 2;
        let mut combos_idx = choose_k_indices(&pair_indices, k);
        if combo_num % 2 == 1 {
            combos_idx.reverse();
        }

        let mut unique_subsets: Vec<Vec<usize>> = Vec::new();
        let mut seen: HashSet<Vec<usize>> = HashSet::new();
        for subset in combos_idx {
            let mut subset_sorted = subset.clone();
            subset_sorted.sort_unstable();
            let mut mirror: Vec<usize> = pair_indices
                .iter()
                .copied()
                .filter(|i| !subset_sorted.contains(i))
                .collect();
            mirror.sort_unstable();

            if seen.contains(&mirror) {
                continue;
            }
            seen.insert(subset_sorted.clone());
            unique_subsets.push(subset_sorted);
        }

        for subset in unique_subsets {
            let mut new_species = combo.species.clone();
            for pidx in subset {
                let a = 2 * pidx;
                let b = a + 1;
                new_species.swap(a, b);
            }
            let mut sp_to_y: HashMap<&str, f64> = HashMap::new();
            for (sp, y) in combo.species.iter().zip(combo.y_raw.iter()) {
                sp_to_y.insert(sp.as_str(), *y);
            }
            let mut new_y = Vec::with_capacity(new_species.len());
            for sp in &new_species {
                new_y.push(*sp_to_y.get(sp.as_str()).unwrap_or(&0.0));
            }

            out.push(ComboJob {
                index: 0,
                combo_label: String::new(),
                combo_tag: String::new(),
                species: new_species,
                y_raw: new_y,
            });
        }
    }

    for (idx, combo) in out.iter_mut().enumerate() {
        let label = format!("combo_{}", idx);
        combo.index = idx;
        combo.combo_label = label.clone();
        combo.combo_tag = make_species_combo_tag(&combo.species, &label);
    }
    Ok(out)
}

fn choose_k_indices(items: &[usize], k: usize) -> Vec<Vec<usize>> {
    fn rec(
        items: &[usize],
        k: usize,
        start: usize,
        cur: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        if start >= items.len() {
            return;
        }
        for i in start..items.len() {
            cur.push(items[i]);
            rec(items, k, i + 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    let mut cur = Vec::new();
    rec(items, k, 0, &mut cur, &mut out);
    out
}

fn checkpoint_paths(output_dir: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf, PathBuf) {
    let dir = output_dir.join("checkpoint");
    let state = dir.join("state.json.gz");
    let cmd = dir.join("command.json");
    let meta = dir.join("meta.txt");
    let runs = dir.join("runs.jsonl");
    (dir, state, cmd, meta, runs)
}

fn normalized_command_value(args: &Args) -> Result<serde_json::Value> {
    let mut v = serde_json::to_value(args)?;
    if let Some(obj) = v.as_object_mut() {
        obj.remove("no_checkpoint");
        obj.remove("force_from_beginning");
        obj.remove("make_sps_plot");
        obj.remove("make_sps_kde_plot");
        obj.remove("make_continuous_plot");
    }
    Ok(v)
}

fn restore_checkpoint_if_available(
    args: &Args,
    output_dir: &Path,
    gene_aggregates: &mut Vec<GeneAggregate>,
    prediction_rows: &mut Vec<PredictionRowOut>,
    total_model_runs: &mut usize,
) -> Result<Option<usize>> {
    let (cp_dir, state_path, cmd_path, meta_path, _runs_path) = checkpoint_paths(output_dir);
    let legacy_uncompressed_state_path = cp_dir.join("state.json");
    let legacy_state_path = cp_dir.join("state_unified.json");
    let legacy_cmd_path = cp_dir.join("command_unified.json");
    if args.force_from_beginning && cp_dir.exists() {
        fs::remove_dir_all(&cp_dir)
            .with_context(|| format!("failed to clear checkpoint dir {}", cp_dir.display()))?;
        return Ok(Some(0));
    }
    let active_state_path = if state_path.exists() {
        state_path.clone()
    } else if legacy_uncompressed_state_path.exists() {
        legacy_uncompressed_state_path
    } else if legacy_state_path.exists() {
        legacy_state_path
    } else {
        return Ok(None);
    };

    let active_cmd_path = if cmd_path.exists() {
        Some(cmd_path.clone())
    } else if legacy_cmd_path.exists() {
        Some(legacy_cmd_path)
    } else {
        None
    };
    if let Some(cmd_to_read) = active_cmd_path {
        let old_cmd: serde_json::Value = serde_json::from_reader(
            File::open(&cmd_to_read)
                .with_context(|| format!("failed to read {}", cmd_to_read.display()))?,
        )?;
        let cur_cmd = normalized_command_value(args)?;
        if old_cmd != cur_cmd {
            bail!(
                "checkpoint exists but command arguments differ; use --force_from_beginning to restart"
            );
        }
    }

    let state_file = File::open(&active_state_path)
        .with_context(|| format!("failed to read {}", active_state_path.display()))?;
    let state: CheckpointState =
        if active_state_path.extension().and_then(|s| s.to_str()) == Some("gz") {
            serde_json::from_reader(GzDecoder::new(state_file)).with_context(|| {
                format!(
                    "failed to decode gzip checkpoint {}",
                    active_state_path.display()
                )
            })?
        } else {
            serde_json::from_reader(state_file).with_context(|| {
                format!("failed to parse checkpoint {}", active_state_path.display())
            })?
        };
    if state.gene_aggregates.len() == gene_aggregates.len() {
        *gene_aggregates = state.gene_aggregates;
    }
    *prediction_rows = state.prediction_rows;
    *total_model_runs = state.total_model_runs;

    // Python parity: meta.txt stores last completed combo index.
    let next_combo_index = if meta_path.exists() {
        match fs::read_to_string(&meta_path)
            .ok()
            .and_then(|s| s.trim().parse::<usize>().ok())
        {
            Some(last_done) => last_done.saturating_add(1),
            None => state.next_combo_index,
        }
    } else {
        state.next_combo_index
    };
    Ok(Some(next_combo_index))
}

fn save_checkpoint_state(
    args: &Args,
    output_dir: &Path,
    next_combo_index: usize,
    total_model_runs: usize,
    gene_aggregates: &[GeneAggregate],
    prediction_rows: &[PredictionRowOut],
    run_audits: &[CheckpointRunAudit],
) -> Result<()> {
    let (cp_dir, state_path, cmd_path, meta_path, runs_path) = checkpoint_paths(output_dir);
    fs::create_dir_all(&cp_dir)
        .with_context(|| format!("failed to create checkpoint dir {}", cp_dir.display()))?;
    let legacy_uncompressed_state_path = cp_dir.join("state.json");
    let legacy_state_path = cp_dir.join("state_unified.json");
    let legacy_cmd_path = cp_dir.join("command_unified.json");

    let state = CheckpointState {
        next_combo_index,
        total_model_runs,
        gene_aggregates: gene_aggregates.to_vec(),
        prediction_rows: prediction_rows.to_vec(),
    };
    {
        let tmp = cp_dir.join("state.json.gz.tmp");
        let tmp_file =
            File::create(&tmp).with_context(|| format!("failed to write {}", tmp.display()))?;
        let mut encoder = GzEncoder::new(tmp_file, Compression::default());
        serde_json::to_writer_pretty(&mut encoder, &state)
            .with_context(|| format!("failed to serialize {}", tmp.display()))?;
        encoder
            .finish()
            .with_context(|| format!("failed to finalize {}", tmp.display()))?;
        fs::rename(&tmp, &state_path).with_context(|| {
            format!(
                "failed to move checkpoint {} -> {}",
                tmp.display(),
                state_path.display()
            )
        })?;
    }

    if !cmd_path.exists() {
        let cmd = normalized_command_value(args)?;
        serde_json::to_writer_pretty(
            File::create(&cmd_path)
                .with_context(|| format!("failed to write {}", cmd_path.display()))?,
            &cmd,
        )?;
    }

    // Python parity: meta.txt stores last completed combo index.
    let last_completed = next_combo_index.saturating_sub(1);
    fs::write(&meta_path, format!("{last_completed}\n"))
        .with_context(|| format!("failed to write {}", meta_path.display()))?;

    if !run_audits.is_empty() {
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&runs_path)
            .with_context(|| format!("failed to open {}", runs_path.display()))?;
        for rec in run_audits {
            serde_json::to_writer(&mut f, rec)?;
            writeln!(f)?;
        }
    }
    if legacy_state_path.exists() {
        let _ = fs::remove_file(&legacy_state_path);
    }
    if legacy_uncompressed_state_path.exists() {
        let _ = fs::remove_file(&legacy_uncompressed_state_path);
    }
    if legacy_cmd_path.exists() {
        let _ = fs::remove_file(&legacy_cmd_path);
    }

    Ok(())
}

fn resolve_combo_alignments<'a>(
    args: &Args,
    combo: &ComboJob,
    default_alignments: &'a [GeneAlignment],
    limited_genes: Option<&HashSet<String>>,
    cache: &'a mut HashMap<PathBuf, Vec<GeneAlignment>>,
) -> Result<&'a [GeneAlignment]> {
    if !args.use_existing_alignments {
        return Ok(default_alignments);
    }
    let base = args
        .canceled_alignments_dir
        .as_ref()
        .ok_or_else(|| anyhow!("--use_existing_alignments requires --canceled_alignments_dir"))?;
    let candidate = base.join(format!("{}-alignments", combo.combo_label));
    let dir = if candidate.is_dir() {
        candidate
    } else {
        base.clone()
    };
    if !dir.is_dir() {
        bail!(
            "existing alignments directory does not exist: {}",
            dir.display()
        );
    }
    let use_dir = if args.auto_convert_to_2line {
        ensure_two_line_dir(&dir, true)?
    } else {
        dir.clone()
    };

    if !cache.contains_key(&use_dir) {
        let loaded = load_alignments(&use_dir, limited_genes, None, None)?;
        cache.insert(use_dir.clone(), loaded);
    }
    Ok(cache
        .get(&use_dir)
        .map(|v| v.as_slice())
        .ok_or_else(|| anyhow!("failed to load combo alignments"))?)
}

fn is_two_line_fasta(path: &Path) -> Result<bool> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut seen_header = false;
    let mut seen_seq_for_header = false;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('>') {
            if seen_header && !seen_seq_for_header {
                bail!(
                    "header without sequence in {} at line {}",
                    path.display(),
                    lineno + 1
                );
            }
            seen_header = true;
            seen_seq_for_header = false;
        } else if seen_header {
            if seen_seq_for_header {
                return Ok(false);
            }
            seen_seq_for_header = true;
        } else if !line.trim().is_empty() {
            bail!(
                "sequence line before first header in {} at line {}",
                path.display(),
                lineno + 1
            );
        }
    }
    if seen_header && !seen_seq_for_header {
        bail!("header without sequence at end of {}", path.display());
    }
    Ok(true)
}

fn ensure_two_line_dir(dir: &Path, recursive: bool) -> Result<PathBuf> {
    let mut bad: Vec<PathBuf> = Vec::new();
    if recursive {
        visit_dirs_collect_fasta(dir, &mut bad, true)?;
    } else {
        for ent in fs::read_dir(dir).with_context(|| format!("unable to read {}", dir.display()))? {
            let p = ent?.path();
            if p.is_file() && is_fasta_path(&p) && !is_two_line_fasta(&p)? {
                bad.push(p);
            }
        }
    }
    if bad.is_empty() {
        return Ok(dir.to_path_buf());
    }

    let parent = dir
        .parent()
        .ok_or_else(|| anyhow!("cannot derive parent for {}", dir.display()))?;
    let base = dir
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("invalid directory name {}", dir.display()))?;
    let target = parent.join(format!("{}_2line", base));
    if target.exists() {
        bail!(
            "cannot auto-convert to 2-line FASTA; target already exists: {}",
            target.display()
        );
    }
    fs::create_dir_all(&target)
        .with_context(|| format!("failed to create {}", target.display()))?;
    convert_dir_to_two_line(dir, &target, recursive)?;
    Ok(target)
}

fn visit_dirs_collect_fasta(dir: &Path, bad: &mut Vec<PathBuf>, recurse: bool) -> Result<()> {
    for ent in fs::read_dir(dir).with_context(|| format!("unable to read {}", dir.display()))? {
        let p = ent?.path();
        if p.is_dir() && recurse {
            visit_dirs_collect_fasta(&p, bad, true)?;
        } else if p.is_file() && is_fasta_path(&p) && !is_two_line_fasta(&p)? {
            bad.push(p);
        }
    }
    Ok(())
}

fn convert_dir_to_two_line(src: &Path, dst: &Path, recursive: bool) -> Result<()> {
    if recursive {
        for ent in fs::read_dir(src).with_context(|| format!("unable to read {}", src.display()))? {
            let p = ent?.path();
            if p.is_dir() {
                let rel = p
                    .strip_prefix(src)
                    .with_context(|| format!("failed to strip prefix {}", p.display()))?;
                let out_sub = dst.join(rel);
                fs::create_dir_all(&out_sub)?;
                convert_dir_to_two_line(&p, &out_sub, true)?;
            } else if p.is_file() && is_fasta_path(&p) {
                convert_single_fasta_to_two_line(
                    &p,
                    &dst.join(
                        p.file_name()
                            .ok_or_else(|| anyhow!("invalid file name {}", p.display()))?,
                    ),
                )?;
            }
        }
    } else {
        for ent in fs::read_dir(src).with_context(|| format!("unable to read {}", src.display()))? {
            let p = ent?.path();
            if p.is_file() && is_fasta_path(&p) {
                convert_single_fasta_to_two_line(
                    &p,
                    &dst.join(
                        p.file_name()
                            .ok_or_else(|| anyhow!("invalid file name {}", p.display()))?,
                    ),
                )?;
            }
        }
    }
    Ok(())
}

fn convert_single_fasta_to_two_line(src: &Path, dst: &Path) -> Result<()> {
    let (seqs, _) = read_fasta_map(src)?;
    let mut out =
        File::create(dst).with_context(|| format!("failed to create {}", dst.display()))?;
    for (species, seq) in seqs {
        writeln!(out, ">{species}")?;
        writeln!(out, "{}", String::from_utf8_lossy(&seq))?;
    }
    Ok(())
}

fn maybe_generate_plots(args: &Args, pred_csv_path: &Path, use_continuous: bool) -> Result<()> {
    if !(args.make_sps_plot || args.make_sps_kde_plot || args.make_continuous_plot) {
        return Ok(());
    }
    if args.make_continuous_plot {
        run_python_plot_command(
            "continuous",
            pred_csv_path,
            &args.output_file_base_name,
            args.min_genes,
            args.pheno_names.clone(),
        )?;
        return Ok(());
    }
    if use_continuous {
        statusln!("Skipping SPS density plots for continuous phenotypes.");
        return Ok(());
    }
    let plot_type = if args.make_sps_plot { "violin" } else { "kde" };
    run_python_plot_command(
        plot_type,
        pred_csv_path,
        &args.output_file_base_name,
        args.min_genes,
        args.pheno_names.clone(),
    )?;
    Ok(())
}

fn run_python_plot_command(
    mode: &str,
    pred_csv_path: &Path,
    title: &str,
    min_genes: usize,
    pheno_names: Option<Vec<String>>,
) -> Result<()> {
    let mut module_args = vec![
        "--mode".to_string(),
        mode.to_string(),
        "--pred_csv".to_string(),
        pred_csv_path.display().to_string(),
        "--title".to_string(),
        title.to_string(),
        "--min_genes".to_string(),
        min_genes.to_string(),
    ];
    if let Some(names) = &pheno_names {
        if names.len() == 2 {
            module_args.push("--pheno_name1".to_string());
            module_args.push(names[0].clone());
            module_args.push("--pheno_name2".to_string());
            module_args.push(names[1].clone());
        }
    }

    if let Ok(override_plotter) = std::env::var("ESL_PSC_PLOTTER") {
        let trimmed = override_plotter.trim();
        if !trimmed.is_empty() {
            let mut cmd = Command::new(trimmed);
            cmd.args(&module_args);
            if let Ok(status) = cmd.status() {
                if status.success() {
                    return Ok(());
                }
                eprintln!(
                    "plot helper '{}' exited with status {}; falling back to python module",
                    trimmed, status
                );
            } else {
                eprintln!(
                    "plot helper '{}' could not be launched; falling back to python module",
                    trimmed
                );
            }
        }
    }

    run_python_module_subcommand("esl_psc_cli.plot_cli", &module_args, &[])
}

fn python_command_candidates() -> Vec<String> {
    let mut candidates: Vec<String> = Vec::new();
    if let Ok(path) = std::env::var("ESL_PSC_PYTHON") {
        let trimmed = path.trim();
        if !trimmed.is_empty() {
            candidates.push(trimmed.to_string());
        }
    }
    if cfg!(windows) {
        candidates.push("python".to_string());
        candidates.push("python3".to_string());
    } else {
        candidates.push("python3".to_string());
        candidates.push("python".to_string());
    }
    dedupe_strings(candidates)
}

fn dedupe_pathbufs(items: Vec<PathBuf>) -> Vec<PathBuf> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = Vec::new();
    for item in items {
        let key = item.to_string_lossy().to_string();
        if seen.insert(key) {
            out.push(item);
        }
    }
    out
}

fn dedupe_strings(items: Vec<String>) -> Vec<String> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = Vec::new();
    for item in items {
        if seen.insert(item.clone()) {
            out.push(item);
        }
    }
    out
}

fn maybe_dump_preprocess_features(
    output_dir: &Path,
    combo_label: &str,
    prep: &PreprocessedData,
) -> Result<()> {
    let enabled = std::env::var("ESL_UNIFIED_DUMP_FEATURES")
        .ok()
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    if !enabled {
        return Ok(());
    }

    let path = output_dir.join(format!("{}_debug_features.txt", combo_label));
    let mut f = File::create(&path)
        .with_context(|| format!("failed to create debug feature dump {}", path.display()))?;
    for fm in &prep.features {
        writeln!(f, "{}", fm.label)?;
    }
    Ok(())
}

fn read_limited_gene_set(path: &Path) -> Result<HashSet<String>> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let stem = Path::new(t)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(t)
            .to_string();
        set.insert(stem);
    }
    Ok(set)
}

fn read_path_file_order(path_file: &Path, alignments_dir: &Path) -> Result<Vec<PathBuf>> {
    let file =
        File::open(path_file).with_context(|| format!("unable to open {}", path_file.display()))?;
    let reader = BufReader::new(file);
    let mut files = Vec::new();

    for line in reader.lines() {
        let raw = line?;
        let t = raw.trim();
        if t.is_empty() {
            continue;
        }

        let mut candidate = PathBuf::from(t);
        if !candidate.is_absolute() {
            candidate = alignments_dir.join(&candidate);
        }
        if !candidate.exists() {
            if candidate.extension().is_none() {
                let with_fas = candidate.with_extension("fas");
                if with_fas.exists() {
                    candidate = with_fas;
                }
            }
        }
        if !candidate.is_file() {
            bail!(
                "alignment listed in --path_file_path was not found: {}",
                candidate.display()
            );
        }
        if !is_fasta_path(&candidate) {
            bail!(
                "non-FASTA entry in --path_file_path: {}",
                candidate.display()
            );
        }
        files.push(candidate);
    }

    if files.is_empty() {
        bail!("--path_file_path contained no alignment entries");
    }

    Ok(files)
}

fn preprocess_dir_name_for_combo(
    args: &Args,
    combo: &ComboJob,
    alignments_dir: &Path,
    is_multimatrix_mode: bool,
) -> String {
    if is_multimatrix_mode {
        return format!("{}_{}", args.output_file_base_name, combo.combo_label);
    }
    if let Some(name) = &args.preprocessed_dir_name {
        return name.clone();
    }
    let align_name = alignments_dir
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("alignments");
    let resp_name = args
        .response_matrix_path
        .as_ref()
        .or(args.response_file.as_ref())
        .and_then(|p| p.file_name().and_then(|s| s.to_str()))
        .unwrap_or("response.txt");
    format!("{}_{}", align_name, resp_name)
}

fn write_preprocessed_data(dir: &Path, prep: &PreprocessedData) -> Result<()> {
    if dir.exists() {
        fs::remove_dir_all(dir)
            .with_context(|| format!("failed clearing preprocess dir {}", dir.display()))?;
    }
    fs::create_dir_all(dir).with_context(|| format!("failed to create {}", dir.display()))?;
    let base = dir
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("invalid preprocess directory name {}", dir.display()))?;

    let feature_path = dir.join(format!("feature_{}.txt", base));
    let mut ff = File::create(&feature_path)
        .with_context(|| format!("failed to create {}", feature_path.display()))?;
    for row in &prep.x {
        let mut first = true;
        for v in row {
            if !first {
                write!(ff, ",")?;
            }
            if (*v - 1.0).abs() < 1e-12 {
                write!(ff, "1.000000")?;
            } else if v.abs() < 1e-12 {
                write!(ff, "0")?;
            } else {
                write!(ff, "{:.6}", v)?;
            }
            first = false;
        }
        writeln!(ff)?;
    }

    let mapping_path = dir.join(format!("feature_mapping_{}.txt", base));
    let mut mf = File::create(&mapping_path)
        .with_context(|| format!("failed to create {}", mapping_path.display()))?;
    writeln!(mf, "0\t")?;
    for (idx, fm) in prep.features.iter().enumerate() {
        writeln!(mf, "{}\t{}", idx + 1, fm.label)?;
    }

    let response_path = dir.join(format!("response_{}.txt", base));
    let mut rf = File::create(&response_path)
        .with_context(|| format!("failed to create {}", response_path.display()))?;
    for y in &prep.y_model {
        writeln!(rf, "{}", format_float_trim(*y))?;
    }

    let stats_path = dir.join(format!("feature_stats_{}.txt", base));
    let mut sf = File::create(&stats_path)
        .with_context(|| format!("failed to create {}", stats_path.display()))?;
    writeln!(sf, "Samples\t{}", prep.species.len())?;
    writeln!(sf, "Features\t{}", prep.features.len())?;

    let field_path = dir.join(format!("field_{}.txt", base));
    let mut fld = File::create(&field_path)
        .with_context(|| format!("failed to create {}", field_path.display()))?;
    for i in 1..=prep.features.len() {
        if i > 1 {
            write!(fld, ",")?;
        }
        write!(fld, "{}", i)?;
    }
    writeln!(fld)?;

    let miss_path = dir.join(format!("missing_seqs_{}.txt", base));
    File::create(&miss_path)
        .with_context(|| format!("failed to create {}", miss_path.display()))?;

    let mut starts = Vec::new();
    let mut ends = Vec::new();
    let mut weights = Vec::new();
    for g in &prep.genes {
        let s = g.feature_start + 1;
        // Keep empty groups so line lengths match legacy preprocess output.
        // For empty ranges, legacy files encode end as start-1 (1-based).
        let e = if g.feature_end >= g.feature_start {
            g.feature_end + 1
        } else {
            g.feature_start
        };
        let len = if e >= s { e - s + 1 } else { 0 };
        starts.push(s.to_string());
        ends.push(e.to_string());
        weights.push(format!("{:.6}", (len as f64).sqrt()));
    }
    let group_path = dir.join(format!("group_indices_{}.txt", base));
    let mut gf = File::create(&group_path)
        .with_context(|| format!("failed to create {}", group_path.display()))?;
    writeln!(gf, "{}", starts.join(","))?;
    writeln!(gf, "{}", ends.join(","))?;
    writeln!(gf, "{}", weights.join(","))?;

    Ok(())
}

fn load_preprocessed_data(
    dir: &Path,
    combo: &ComboJob,
    combo_alignments: &[GeneAlignment],
    continuous: bool,
) -> Result<PreprocessedData> {
    if !dir.is_dir() {
        bail!(
            "--use_existing_preprocess requested but directory not found: {}",
            dir.display()
        );
    }
    let base = dir
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("invalid preprocess directory name {}", dir.display()))?;
    let feature_path = dir.join(format!("feature_{}.txt", base));
    let mapping_path = dir.join(format!("feature_mapping_{}.txt", base));
    let response_path = dir.join(format!("response_{}.txt", base));

    let mut features: Vec<FeatureMeta> = Vec::new();
    for line in BufReader::new(
        File::open(&mapping_path)
            .with_context(|| format!("failed to open {}", mapping_path.display()))?,
    )
    .lines()
    {
        let line = line?;
        let mut parts = line.splitn(2, '\t');
        let idx_part = parts.next().unwrap_or("").trim();
        let label = parts.next().unwrap_or("").trim();
        if idx_part.is_empty() || label.is_empty() || !idx_part.chars().all(|c| c.is_ascii_digit())
        {
            continue;
        }
        if idx_part == "0" {
            continue;
        }
        let (gene_name, pos, aa) = parse_feature_label(label)?;
        features.push(FeatureMeta {
            label: label.to_string(),
            gene_idx: 0,
            position: pos,
            aa,
        });
        let _ = gene_name;
    }

    let mut x: Vec<Vec<f64>> = Vec::new();
    for line in BufReader::new(
        File::open(&feature_path)
            .with_context(|| format!("failed to open {}", feature_path.display()))?,
    )
    .lines()
    {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let mut row = Vec::new();
        for tok in t.split(',') {
            let v = tok.trim().parse::<f64>().with_context(|| {
                format!(
                    "invalid numeric token '{}' in {}",
                    tok,
                    feature_path.display()
                )
            })?;
            row.push(v);
        }
        x.push(row);
    }
    if x.is_empty() {
        bail!("no samples parsed from {}", feature_path.display());
    }
    if !features.is_empty() && x[0].len() != features.len() {
        bail!(
            "feature column count {} does not match mapping count {} in {}",
            x[0].len(),
            features.len(),
            dir.display()
        );
    }

    let mut y_file: Vec<f64> = Vec::new();
    for line in BufReader::new(
        File::open(&response_path)
            .with_context(|| format!("failed to open {}", response_path.display()))?,
    )
    .lines()
    {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        y_file.push(
            t.parse::<f64>().with_context(|| {
                format!("invalid response '{}' in {}", t, response_path.display())
            })?,
        );
    }

    let y_model = if continuous {
        if y_file.len() == combo.species.len() {
            y_file
        } else {
            combo.y_raw.clone()
        }
    } else {
        binary_pair_labels(combo.species.len())
    };

    let var_counts = compute_var_site_counts(combo_alignments, &combo.species);
    let mut feature_ranges: HashMap<String, (usize, usize)> = HashMap::new();
    if !features.is_empty() {
        let mut cur_gene = parse_feature_label(&features[0].label)?.0;
        let mut cur_start = 0usize;
        for j in 0..features.len() {
            let (gene_name, _p, _a) = parse_feature_label(&features[j].label)?;
            if gene_name != cur_gene {
                feature_ranges.insert(cur_gene.clone(), (cur_start, j - 1));
                cur_gene = gene_name;
                cur_start = j;
            }
        }
        feature_ranges.insert(cur_gene, (cur_start, features.len() - 1));
    }

    let mut genes: Vec<GeneMeta> = Vec::with_capacity(combo_alignments.len());
    for (gidx, g) in combo_alignments.iter().enumerate() {
        let (start, end) = if let Some((s, e)) = feature_ranges.get(&g.name) {
            for k in *s..=*e {
                features[k].gene_idx = gidx;
            }
            (*s, *e)
        } else {
            // Represent empty feature range as end < start so solver paths skip it.
            (1, 0)
        };
        genes.push(GeneMeta {
            name: g.name.clone(),
            feature_start: start,
            feature_end: end,
            var_site_count: *var_counts.get(&g.name).unwrap_or(&0),
        });
    }

    Ok(PreprocessedData {
        species: combo.species.clone(),
        y_model,
        x,
        features,
        genes,
    })
}

fn parse_feature_label(label: &str) -> Result<(String, usize, u8)> {
    let mut parts = label.rsplitn(3, '_');
    let aa = parts
        .next()
        .ok_or_else(|| anyhow!("invalid feature label '{}'", label))?;
    let pos = parts
        .next()
        .ok_or_else(|| anyhow!("invalid feature label '{}'", label))?;
    let gene = parts
        .next()
        .ok_or_else(|| anyhow!("invalid feature label '{}'", label))?;
    let pos = pos
        .parse::<usize>()
        .with_context(|| format!("invalid feature position in '{}'", label))?;
    let aa = aa
        .as_bytes()
        .first()
        .copied()
        .ok_or_else(|| anyhow!("invalid feature residue in '{}'", label))?;
    Ok((gene.to_string(), pos, aa))
}

fn compute_var_site_counts(
    alignments: &[GeneAlignment],
    combo_species: &[String],
) -> HashMap<String, usize> {
    let mut out = HashMap::new();
    for gene in alignments {
        let mut seqs = Vec::with_capacity(combo_species.len());
        for sp in combo_species {
            if let Some(s) = gene.seqs.get(sp) {
                seqs.push(s.clone());
            } else {
                seqs.push(vec![b'-'; gene.seq_len]);
            }
        }
        out.insert(gene.name.clone(), count_var_sites_python(&seqs));
    }
    out
}

fn load_alignments(
    alignments_dir: &Path,
    limited_genes: Option<&HashSet<String>>,
    ordered_fasta_files: Option<&[PathBuf]>,
    progress_label: Option<&str>,
) -> Result<Vec<GeneAlignment>> {
    let fasta_files: Vec<PathBuf> = if let Some(ordered) = ordered_fasta_files {
        ordered.to_vec()
    } else {
        // Preserve directory iteration order to match Python's os.listdir()
        // behavior used to build preprocess paths.txt in the legacy pipeline.
        fs::read_dir(alignments_dir)
            .with_context(|| format!("unable to read {}", alignments_dir.display()))?
            .filter_map(|e| e.ok().map(|x| x.path()))
            .filter(|p| p.is_file() && is_fasta_path(p))
            .collect()
    };
    let fasta_files: Vec<PathBuf> = if let Some(allowed) = limited_genes {
        fasta_files
            .into_iter()
            .filter(|fasta| {
                fasta.file_stem()
                    .and_then(|s| s.to_str())
                    .map(|gene_name| allowed.contains(gene_name))
                    .unwrap_or(false)
            })
            .collect()
    } else {
        fasta_files
    };

    let total_files = fasta_files.len();
    if let Some(label) = progress_label {
        statusln!("Loading {} alignment files: {}", label, total_files);
    }
    let progress_interval = ((total_files + 49) / 50).max(1);

    let mut out = Vec::with_capacity(total_files);
    for (idx, fasta) in fasta_files.into_iter().enumerate() {
        let gene_name = fasta
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("invalid FASTA file name: {}", fasta.display()))?
            .to_string();

        let (seq_map, seq_len) = read_fasta_map(&fasta)?;
        out.push(GeneAlignment {
            name: gene_name,
            seq_len,
            seqs: seq_map,
        });
        if let Some(label) = progress_label {
            let current = idx + 1;
            if current == 1 || current == total_files || current % progress_interval == 0 {
                statusln!("Loaded {} alignment file {} of {}", label, current, total_files);
            }
        }
    }

    Ok(out)
}

fn read_fasta_map(path: &Path) -> Result<(HashMap<String, Vec<u8>>, usize)> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut seq_map: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_id = String::new();
    let mut current_seq: Vec<u8> = Vec::new();
    let mut seq_len: Option<usize> = None;

    for line in reader.lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(rest) = t.strip_prefix('>') {
            if !current_id.is_empty() {
                let current_len = current_seq.len();
                if let Some(expected) = seq_len {
                    if current_len != expected {
                        bail!(
                            "sequence length mismatch in {} for species {}",
                            path.display(),
                            current_id
                        );
                    }
                } else {
                    seq_len = Some(current_len);
                }
                seq_map.insert(
                    std::mem::take(&mut current_id),
                    std::mem::take(&mut current_seq),
                );
            }
            current_id = rest.trim().to_string();
        } else if !current_id.is_empty() {
            for b in t.as_bytes() {
                let c = if *b == b'?' {
                    b'-'
                } else {
                    b.to_ascii_uppercase()
                };
                current_seq.push(c);
            }
        }
    }

    if !current_id.is_empty() {
        let current_len = current_seq.len();
        if let Some(expected) = seq_len {
            if current_len != expected {
                bail!(
                    "sequence length mismatch in {} for species {}",
                    path.display(),
                    current_id
                );
            }
        } else {
            seq_len = Some(current_len);
        }
        seq_map.insert(current_id, current_seq);
    }

    if seq_map.is_empty() {
        bail!("no sequences parsed from {}", path.display());
    }

    Ok((seq_map, seq_len.unwrap_or(0)))
}

fn write_two_line_fasta(path: &Path, species: &[String], seqs: &[Vec<u8>]) -> Result<()> {
    if species.len() != seqs.len() {
        bail!(
            "cannot write FASTA {}; species/sequence count mismatch",
            path.display()
        );
    }
    let mut out =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    for (sp, seq) in species.iter().zip(seqs.iter()) {
        writeln!(out, ">{sp}")?;
        writeln!(out, "{}", String::from_utf8_lossy(seq))?;
    }
    Ok(())
}

fn read_species_phenotypes(path: &Path) -> Result<PhenotypeInfo> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut values = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let parts: Vec<&str> = t.split(',').map(|s| s.trim()).collect();
        if parts.len() != 2 || parts[0].is_empty() {
            continue;
        }
        if parts[0].eq_ignore_ascii_case("species") {
            continue;
        }
        let val = match parts[1].parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        values.insert(parts[0].to_string(), val);
    }

    let is_binary = !values.is_empty() && values.values().all(|v| is_binary_pheno_value(*v));

    Ok(PhenotypeInfo { values, is_binary })
}

fn is_binary_pheno_value(v: f64) -> bool {
    (v - 1.0).abs() < 1e-12 || (v + 1.0).abs() < 1e-12 || v.abs() < 1e-12
}

fn binary_pair_labels(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect()
}

fn build_combo_jobs(args: &Args, phenotype_info: Option<&PhenotypeInfo>) -> Result<Vec<ComboJob>> {
    if let Some(path) = args
        .response_matrix_path
        .as_ref()
        .or(args.response_file.as_ref())
    {
        let (species, y_raw) = read_response_matrix(path)?;
        let label = "combo_0".to_string();
        return Ok(vec![ComboJob {
            index: 0,
            combo_label: label.clone(),
            combo_tag: make_species_combo_tag(&species, &label),
            species,
            y_raw,
        }]);
    }

    if let Some(dir) = &args.response_dir {
        let mut files: Vec<PathBuf> = fs::read_dir(dir)
            .with_context(|| format!("unable to read response_dir {}", dir.display()))?
            .filter_map(|e| e.ok().map(|x| x.path()))
            .filter(|p| {
                p.is_file()
                    && p.extension()
                        .and_then(|s| s.to_str())
                        .map(|s| s.eq_ignore_ascii_case("txt"))
                        .unwrap_or(false)
            })
            .collect();
        files.sort();

        if files.is_empty() {
            bail!("no .txt response matrices found in {}", dir.display());
        }

        let mut jobs = Vec::new();
        for (idx, file) in files.iter().enumerate() {
            let (species, y_raw) = read_response_matrix(file)?;
            let label = format!("combo_{}", idx);
            jobs.push(ComboJob {
                index: idx,
                combo_label: label.clone(),
                combo_tag: make_species_combo_tag(&species, &label),
                species,
                y_raw,
            });
        }
        return Ok(jobs);
    }

    let groups_file = args
        .species_groups_file
        .as_ref()
        .ok_or_else(|| anyhow!("species groups file required"))?;

    let groups = parse_species_groups(groups_file)?;
    let combos = expand_species_groups(&groups);

    let mut jobs = Vec::new();
    for (idx, species_combo) in combos.into_iter().enumerate() {
        let y_raw = if args.use_continuous_phenotypes {
            let pheno = phenotype_info.ok_or_else(|| {
                anyhow!(
                    "--use_continuous_phenotypes with --species_groups_file requires --species_pheno_path"
                )
            })?;
            let mut vals = Vec::with_capacity(species_combo.len());
            for sp in &species_combo {
                let v = pheno
                    .values
                    .get(sp)
                    .ok_or_else(|| anyhow!("species '{}' missing from phenotype file", sp))?;
                vals.push(*v);
            }
            vals
        } else {
            species_combo
                .iter()
                .enumerate()
                .map(|(i, _)| if i % 2 == 0 { 1.0 } else { -1.0 })
                .collect()
        };

        let label = format!("combo_{}", idx);
        jobs.push(ComboJob {
            index: idx,
            combo_label: label.clone(),
            combo_tag: make_species_combo_tag(&species_combo, &label),
            species: species_combo,
            y_raw,
        });
    }

    Ok(jobs)
}

fn parse_species_groups(path: &Path) -> Result<Vec<Vec<String>>> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut groups = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let members: Vec<String> = t.split(',').map(|s| s.trim().to_string()).collect();
        if members.iter().any(|s| s.is_empty()) {
            bail!(
                "invalid species groups format in {} at line {}",
                path.display(),
                i + 1
            );
        }
        groups.push(members);
    }

    if groups.is_empty() {
        bail!("species groups file is empty: {}", path.display());
    }
    if groups.len() % 2 != 0 {
        bail!(
            "species groups file must have an even number of lines: {}",
            path.display()
        );
    }

    Ok(groups)
}

fn expand_species_groups(groups: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut combos: Vec<Vec<String>> = vec![Vec::new()];
    for group in groups {
        let mut next = Vec::new();
        for combo in &combos {
            for member in group {
                let mut c = combo.clone();
                c.push(member.clone());
                next.push(c);
            }
        }
        combos = next;
    }
    combos
}

fn make_species_combo_tag(species: &[String], label: &str) -> String {
    if species.len() > 12 {
        return label.to_string();
    }
    let use_underscore = species.first().map(|s| s.contains('_')).unwrap_or(false);

    let parts: Vec<String> = species
        .iter()
        .map(|s| {
            if use_underscore {
                s.split('_').nth(1).unwrap_or(s).chars().take(3).collect()
            } else {
                s.chars().take(3).collect()
            }
        })
        .collect();

    parts.join(".")
}

fn validate_combo_species_against_alignments(
    combos: &[ComboJob],
    alignments: &[GeneAlignment],
) -> Result<()> {
    let mut present = HashSet::new();
    for gene in alignments {
        for sp in gene.seqs.keys() {
            present.insert(sp.clone());
        }
    }

    let mut missing = BTreeSet::new();
    for combo in combos {
        for sp in &combo.species {
            if !present.contains(sp) {
                missing.insert(sp.clone());
            }
        }
    }

    if !missing.is_empty() {
        bail!(
            "species from response/groups were not found in alignments: {}",
            missing.into_iter().collect::<Vec<String>>().join(", ")
        );
    }

    Ok(())
}

fn read_response_matrix(path: &Path) -> Result<(Vec<String>, Vec<f64>)> {
    let file = File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut species = Vec::new();
    let mut values = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }

        let mut parts: Vec<&str> = if t.contains('\t') {
            t.split('\t').map(|s| s.trim()).collect()
        } else if t.contains(',') {
            t.split(',').map(|s| s.trim()).collect()
        } else {
            t.split_whitespace().collect()
        };

        if parts.len() < 2 {
            continue;
        }

        let name = parts.remove(0).to_string();
        let val = parts
            .remove(0)
            .parse::<f64>()
            .with_context(|| format!("invalid phenotype value for species {}", name))?;
        if val.abs() < 1e-12 {
            continue;
        }
        species.push(name);
        values.push(val);
    }

    if species.is_empty() {
        bail!(
            "no species/phenotype entries were parsed from {}",
            path.display()
        );
    }

    Ok((species, values))
}

fn preprocess_combo_alignments(
    alignments: &[GeneAlignment],
    combo: &ComboJob,
    continuous: bool,
    apply_gap_cancel: bool,
    randomize_pairs: bool,
    preserve_alignments_dir: Option<&Path>,
    args: &Args,
) -> Result<PreprocessedData> {
    let mut x: Vec<Vec<f64>> = vec![Vec::new(); combo.species.len()];
    let mut features: Vec<FeatureMeta> = Vec::new();
    let mut genes: Vec<GeneMeta> = Vec::with_capacity(alignments.len());

    // Process genes in parallel chunks while merging each chunk in input order,
    // preserving deterministic feature/gene ordering.
    let chunk_size = rayon::current_num_threads().max(1);
    for (chunk_idx, gene_chunk) in alignments.chunks(chunk_size).enumerate() {
        let chunk_start = chunk_idx * chunk_size;
        let chunk_results: Vec<GenePreprocessResult> = if gene_chunk.len() == 1 {
            vec![preprocess_one_gene(
                &gene_chunk[0],
                combo,
                apply_gap_cancel,
                randomize_pairs,
                preserve_alignments_dir,
                args,
            )?]
        } else {
            gene_chunk
                .par_iter()
                .map(|gene| {
                    preprocess_one_gene(
                        gene,
                        combo,
                        apply_gap_cancel,
                        randomize_pairs,
                        preserve_alignments_dir,
                        args,
                    )
                })
                .collect::<Result<Vec<_>>>()?
        };

        for (offset, gene_out) in chunk_results.into_iter().enumerate() {
            let gene_idx = chunk_start + offset;
            let feature_start = features.len();
            let gene_name = gene_out.gene_name;
            for col in gene_out.columns {
                for (row, v) in x.iter_mut().zip(col.values.into_iter()) {
                    row.push(v);
                }
                let label = format!("{}_{}_{}", gene_name, col.position, col.aa as char);
                features.push(FeatureMeta {
                    label,
                    gene_idx,
                    position: col.position,
                    aa: col.aa,
                });
            }
            let feature_end = features.len().saturating_sub(1);
            genes.push(GeneMeta {
                name: gene_name,
                feature_start,
                feature_end,
                var_site_count: gene_out.var_site_count,
            });
        }
    }

    let y_model = if continuous {
        combo.y_raw.clone()
    } else {
        binary_pair_labels(combo.species.len())
    };

    Ok(PreprocessedData {
        species: combo.species.clone(),
        y_model,
        x,
        features,
        genes,
    })
}

fn preprocess_one_gene(
    gene: &GeneAlignment,
    combo: &ComboJob,
    apply_gap_cancel: bool,
    randomize_pairs: bool,
    preserve_alignments_dir: Option<&Path>,
    args: &Args,
) -> Result<GenePreprocessResult> {
    let mut seqs = Vec::with_capacity(combo.species.len());
    let mut missing_flags = Vec::with_capacity(combo.species.len());

    for sp in &combo.species {
        if let Some(seq) = gene.seqs.get(sp) {
            seqs.push(seq.clone());
            missing_flags.push(false);
        } else {
            seqs.push(vec![b'-'; gene.seq_len]);
            missing_flags.push(true);
        }
    }

    if args.nix_full_deletions && missing_flags.iter().any(|v| *v) {
        return Ok(GenePreprocessResult {
            gene_name: gene.name.clone(),
            var_site_count: 0,
            columns: Vec::new(),
        });
    }

    if apply_gap_cancel {
        let outgroup = args
            .outgroup_species
            .as_ref()
            .and_then(|sp| gene.seqs.get(sp).cloned());
        apply_gap_cancellation(
            &mut seqs,
            &missing_flags,
            outgroup.as_ref(),
            args.min_pairs,
            args.cancel_only_partner,
            args.cancel_tri_allelic,
        );
    }
    if randomize_pairs {
        apply_pair_randomization(&mut seqs);
    }
    if let Some(dir) = preserve_alignments_dir {
        let out_path = dir.join(format!("{}.fas", gene.name));
        write_two_line_fasta(&out_path, &combo.species, &seqs)?;
    }

    let (var_site_count, columns) = build_gene_feature_columns(&seqs);
    Ok(GenePreprocessResult {
        gene_name: gene.name.clone(),
        var_site_count,
        columns,
    })
}

fn build_gene_feature_columns(seqs: &[Vec<u8>]) -> (usize, Vec<GeneFeatureColumn>) {
    if seqs.is_empty() {
        return (0, Vec::new());
    }
    let seq_len = seqs[0].len();
    let n_species = seqs.len();
    let mut var_site_count = 0usize;
    let mut columns: Vec<GeneFeatureColumn> = Vec::new();

    for pos in 0..seq_len {
        let mut counts = [0usize; 256];
        let mut seen = [false; 256];
        let mut residues: Vec<u8> = Vec::with_capacity(8);

        for seq in seqs {
            let b = seq[pos];
            let idx = b as usize;
            counts[idx] += 1;
            if !seen[idx] {
                seen[idx] = true;
                residues.push(b);
            }
        }
        residues.sort_unstable();

        let mut non_gap_total = 0usize;
        let mut max_count = 0usize;
        let mut unique_non_gap = 0usize;
        for &b in &residues {
            if b == b'-' {
                continue;
            }
            let c = counts[b as usize];
            unique_non_gap += 1;
            non_gap_total += c;
            if c > max_count {
                max_count = c;
            }
        }

        if unique_non_gap < 2 {
            continue;
        }
        if unique_non_gap == non_gap_total {
            continue;
        }
        if max_count == non_gap_total.saturating_sub(1) {
            continue;
        }
        var_site_count += 1;

        // Match preprocess ignore-singleton behavior (countThreshold=2).
        if non_gap_total.saturating_sub(max_count) < 2 {
            continue;
        }

        for &base in &residues {
            if base == b'-' {
                continue;
            }
            let feature_sum = counts[base as usize];
            if feature_sum < 2 {
                continue;
            }
            let mut values = Vec::with_capacity(n_species);
            for seq in seqs {
                values.push(if seq[pos] == base { 1.0 } else { 0.0 });
            }
            columns.push(GeneFeatureColumn {
                position: pos,
                aa: base,
                values,
            });
        }
    }

    (var_site_count, columns)
}

fn apply_pair_randomization(seqs: &mut [Vec<u8>]) {
    if seqs.len() < 2 || seqs.len() % 2 != 0 {
        return;
    }
    let seq_len = seqs[0].len();
    let mut rng = rand::thread_rng();
    for pos in 0..seq_len {
        for pair_idx in 0..(seqs.len() / 2) {
            if rng.gen_bool(0.5) {
                let a = 2 * pair_idx;
                let b = a + 1;
                let tmp = seqs[a][pos];
                seqs[a][pos] = seqs[b][pos];
                seqs[b][pos] = tmp;
            }
        }
    }
}

fn apply_gap_cancellation(
    seqs: &mut [Vec<u8>],
    missing_flags: &[bool],
    outgroup: Option<&Vec<u8>>,
    min_pairs: usize,
    cancel_only_partner: bool,
    cancel_tri_allelic: bool,
) {
    if seqs.is_empty() {
        return;
    }
    let seq_len = seqs[0].len();

    // Cancel everything if too few intact input pairs are present for this gene.
    if seqs.len() % 2 == 0 {
        let mut intact_pairs = 0usize;
        for pair in missing_flags.chunks(2) {
            if pair.len() == 2 && !pair[0] && !pair[1] {
                intact_pairs += 1;
            }
        }
        if intact_pairs < min_pairs {
            for seq in seqs.iter_mut() {
                for c in seq.iter_mut() {
                    *c = b'-';
                }
            }
            return;
        }
    }

    for pos in 0..seq_len {
        let has_gap = seqs.iter().any(|s| s[pos] == b'-');

        if cancel_only_partner && has_gap {
            let mut pairs_left = seqs.len() / 2;
            for pair_idx in 0..(seqs.len() / 2) {
                let a = seqs[2 * pair_idx][pos];
                let b = seqs[2 * pair_idx + 1][pos];
                if a == b'-' || b == b'-' {
                    seqs[2 * pair_idx][pos] = b'-';
                    seqs[2 * pair_idx + 1][pos] = b'-';
                    pairs_left = pairs_left.saturating_sub(1);
                }
            }
            if pairs_left < min_pairs {
                for seq in seqs.iter_mut() {
                    seq[pos] = b'-';
                }
                continue;
            }
        } else if has_gap {
            for seq in seqs.iter_mut() {
                seq[pos] = b'-';
            }
            continue;
        }

        if let Some(out) = outgroup {
            if pos < out.len() {
                let out_res = out[pos];
                if out_res != b'-' {
                    let mut mismatch = false;
                    for ctrl_idx in (1..seqs.len()).step_by(2) {
                        if seqs[ctrl_idx][pos] != out_res {
                            mismatch = true;
                            break;
                        }
                    }
                    if mismatch {
                        for seq in seqs.iter_mut() {
                            seq[pos] = b'-';
                        }
                        continue;
                    }
                }
            }
        }

        if cancel_tri_allelic && seqs.len() == 4 {
            let mut uniq = BTreeSet::new();
            for seq in seqs.iter() {
                uniq.insert(seq[pos]);
            }
            if uniq.len() == 3 {
                for seq in seqs.iter_mut() {
                    seq[pos] = b'-';
                }
            }
        }
    }
}

fn count_var_sites_python(seqs: &[Vec<u8>]) -> usize {
    if seqs.is_empty() {
        return 0;
    }
    let len = seqs[0].len();
    let mut total = 0usize;

    for i in 0..len {
        let mut counts: HashMap<u8, usize> = HashMap::new();
        for s in seqs {
            let c = s[i];
            if c == b'-' {
                continue;
            }
            *counts.entry(c).or_insert(0) += 1;
        }
        if counts.is_empty() || counts.len() == 1 {
            continue;
        }
        let sum_counts: usize = counts.values().sum();
        if counts.len() == sum_counts {
            continue;
        }
        let max_count = counts.values().copied().max().unwrap_or(0);
        if max_count == sum_counts.saturating_sub(1) {
            continue;
        }
        total += 1;
    }

    total
}

fn build_prediction_design(
    prediction_alignments: &[GeneAlignment],
    features: &[FeatureMeta],
    genes: &[GeneMeta],
    input_species: &[String],
    phenotype_info: Option<&PhenotypeInfo>,
) -> Result<PredictionDesign> {
    let input_set: HashSet<&str> = input_species.iter().map(String::as_str).collect();

    let mut all_species: BTreeSet<String> = BTreeSet::new();
    for gene in prediction_alignments {
        for sp in gene.seqs.keys() {
            all_species.insert(sp.clone());
        }
    }

    let mut species = Vec::new();
    let mut true_values = Vec::new();
    for sp in all_species {
        if input_set.contains(sp.as_str()) {
            continue;
        }
        if let Some(ph) = phenotype_info {
            if !ph.values.contains_key(&sp) {
                continue;
            }
            true_values.push(ph.values.get(&sp).copied());
        } else {
            true_values.push(None);
        }
        species.push(sp);
    }

    let gene_lookup: HashMap<&str, &GeneAlignment> = prediction_alignments
        .iter()
        .map(|g| (g.name.as_str(), g))
        .collect();

    let mut feature_hit_rows = vec![Vec::new(); features.len()];
    for (row_idx, sp) in species.iter().enumerate() {
        for (j, fm) in features.iter().enumerate() {
            let gene_name = genes[fm.gene_idx].name.as_str();
            if let Some(gene) = gene_lookup.get(gene_name) {
                if let Some(seq) = gene.seqs.get(sp) {
                    if fm.position < seq.len() && seq[fm.position] == fm.aa {
                        feature_hit_rows[j].push(row_idx);
                    }
                }
            }
        }
    }

    Ok(PredictionDesign {
        species,
        feature_hit_rows,
        true_values,
    })
}

fn build_feature_hit_rows_from_dense(x: &[Vec<f64>]) -> Vec<Vec<usize>> {
    if x.is_empty() {
        return Vec::new();
    }
    let p = x[0].len();
    let mut rows = vec![Vec::new(); p];
    for (i, row) in x.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            if *v != 0.0 {
                rows[j].push(i);
            }
        }
    }
    rows
}

fn append_prediction_rows(
    out_rows: &mut Vec<PredictionRowOut>,
    prep: &PreprocessedData,
    input_feature_hit_rows: &[Vec<usize>],
    pred: &PredictionDesign,
    run: &ModelResult,
    combo_tag: &str,
    phenotype_info: Option<&PhenotypeInfo>,
) {
    if pred.species.is_empty() {
        return;
    }

    // Python parity: species predictions are emitted only for species that
    // matched at least one selected non-zero feature in that run.
    // Intercept and RMSE are then applied only to species present in the score map.
    let (pred_scores, input_rmse) =
        python_style_prediction_scores(prep, input_feature_hit_rows, pred, run);

    for (i, score) in pred_scores {
        let true_phenotype = if let Some(ph) = phenotype_info {
            match pred.true_values[i] {
                Some(v) if ph.is_binary && v.abs() < 1e-12 => Some(String::new()),
                Some(v) => Some(format_float_trim(v)),
                None => Some(String::new()),
            }
        } else {
            None
        };

        out_rows.push(PredictionRowOut {
            species_combo: combo_tag.to_string(),
            lambda1: run.lambda1,
            lambda2: run.lambda2,
            penalty_term: run.penalty_term,
            num_genes: run.num_genes,
            input_rmse,
            species: pred.species[i].clone(),
            sps: score,
            true_phenotype,
        });
    }
}

fn python_style_prediction_scores(
    prep: &PreprocessedData,
    input_feature_hit_rows: &[Vec<usize>],
    pred: &PredictionDesign,
    run: &ModelResult,
) -> (Vec<(usize, f64)>, f64) {
    let active_features: Vec<usize> = run
        .beta
        .iter()
        .enumerate()
        .filter_map(|(j, w)| if *w != 0.0 { Some(j) } else { None })
        .collect();

    let mut input_scores = vec![0.0_f64; prep.species.len()];
    let mut input_touched = vec![false; prep.species.len()];
    let mut pred_scores_dense = vec![0.0_f64; pred.species.len()];
    let mut pred_touched_mask = vec![false; pred.species.len()];

    for &j in &active_features {
        let w = run.beta[j];

        if let Some(rows) = input_feature_hit_rows.get(j) {
            for &row in rows {
                input_scores[row] += w;
                input_touched[row] = true;
            }
        }
        if let Some(rows) = pred.feature_hit_rows.get(j) {
            for &row in rows {
                pred_scores_dense[row] += w;
                pred_touched_mask[row] = true;
            }
        }
    }

    // Python adds intercept only to species keys already present in the map.
    for i in 0..input_scores.len() {
        if input_touched[i] {
            input_scores[i] += run.intercept;
        }
    }
    for i in 0..pred_scores_dense.len() {
        if pred_touched_mask[i] {
            pred_scores_dense[i] += run.intercept;
        }
    }

    let mut pred_touched: Vec<usize> = Vec::new();
    for i in 0..pred.species.len() {
        if pred_touched_mask[i] {
            pred_touched.push(i);
        }
    }

    // Python computes RMSE over all input species, defaulting missing scores to 0.
    let mut sum_sq = 0.0_f64;
    let n = prep.species.len().max(1);
    for (i, _) in prep.species.iter().enumerate() {
        let observed = if input_touched[i] {
            input_scores[i]
        } else {
            0.0
        };
        let expected = prep.y_model[i];
        sum_sq += (expected - observed).powi(2);
    }
    let input_rmse = (sum_sq / (n as f64)).sqrt();

    let pred_scores = pred_touched
        .into_iter()
        .map(|i| (i, pred_scores_dense[i]))
        .collect();

    (pred_scores, input_rmse)
}

fn update_gene_stats_for_run(
    run: &ModelResult,
    combo_highest_gss: &mut [f64],
    combo_best_rank: &mut [Option<usize>],
    gene_aggregates: &mut [GeneAggregate],
) {
    let mut ranked: Vec<(usize, f64)> = run
        .gene_gss
        .iter()
        .enumerate()
        .filter(|(_, v)| **v > 0.0)
        .map(|(idx, v)| (idx, *v))
        .collect();
    // Python uses a stable sort over per-gene GSS values; equal scores preserve
    // first-seen order. Reproduce that behavior by using gene index as a
    // deterministic tie-breaker.
    ranked.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.0.cmp(&b.0))
    });

    for (gidx, gss) in &ranked {
        if *gss > combo_highest_gss[*gidx] {
            combo_highest_gss[*gidx] = *gss;
        }
        if *gss > gene_aggregates[*gidx].single_highest_gss {
            gene_aggregates[*gidx].single_highest_gss = *gss;
        }
    }

    for (rank_idx, (gidx, _)) in ranked.iter().enumerate() {
        let rank = rank_idx + 1;
        update_best_rank(&mut combo_best_rank[*gidx], rank);
        update_best_rank(&mut gene_aggregates[*gidx].single_best_rank, rank);
    }

    for ss in &run.selected_sites {
        let entry = gene_aggregates[ss.gene_idx]
            .selected_sites
            .entry(ss.position)
            .or_insert(ss.pss);
        if ss.pss > *entry {
            *entry = ss.pss;
        }
    }
}

fn finalize_combo_multimatrix_stats(
    combo_highest_gss: &[f64],
    combo_best_rank: &[Option<usize>],
    top_rank_threshold: f64,
    gene_aggregates: &mut [GeneAggregate],
) {
    for (idx, agg) in gene_aggregates.iter_mut().enumerate() {
        if let Some(rank) = combo_best_rank[idx] {
            agg.num_combos_ranked += 1;
            if rank as f64 <= top_rank_threshold {
                agg.num_combos_ranked_top += 1;
            }
            if combo_highest_gss[idx] > agg.highest_ever_gss {
                agg.highest_ever_gss = combo_highest_gss[idx];
            }
            update_best_rank(&mut agg.best_ever_rank, rank);
        }
    }
}

fn update_best_rank(slot: &mut Option<usize>, rank: usize) {
    match *slot {
        Some(old) if rank >= old => {}
        _ => *slot = Some(rank),
    }
}

fn intercept_only_model(
    prep: &PreprocessedData,
    lambda1: f64,
    lambda2: f64,
    penalty: f64,
    continuous: bool,
) -> ModelResult {
    let intercept = initial_intercept(&prep.y_model, continuous);
    let scores = vec![intercept; prep.species.len()];
    let rmse = compute_rmse(&scores, &prep.y_model);

    ModelResult {
        lambda1,
        lambda2,
        penalty_term: penalty,
        beta: Vec::new(),
        intercept,
        rmse,
        num_genes: 0,
        gene_gss: vec![0.0; prep.genes.len()],
        selected_sites: Vec::new(),
    }
}

fn write_predictions_csv(
    output_dir: &Path,
    base_name: &str,
    rows: &[PredictionRowOut],
    include_true_phenotype: bool,
) -> Result<PathBuf> {
    let path = output_dir.join(format!("{}_species_predictions.csv", base_name));
    let mut wtr =
        Writer::from_path(&path).with_context(|| format!("failed to create {}", path.display()))?;

    if include_true_phenotype {
        wtr.write_record([
            "species_combo",
            "lambda1",
            "lambda2",
            "penalty_term",
            "num_genes",
            "input_RMSE",
            "species",
            "SPS",
            "true_phenotype",
        ])?;
    } else {
        wtr.write_record([
            "species_combo",
            "lambda1",
            "lambda2",
            "penalty_term",
            "num_genes",
            "input_RMSE",
            "species",
            "SPS",
        ])?;
    }

    for row in rows {
        if include_true_phenotype {
            wtr.write_record([
                row.species_combo.clone(),
                format_float_trim(row.lambda1),
                format_float_trim(row.lambda2),
                format_float_trim(row.penalty_term),
                row.num_genes.to_string(),
                format_float_trim(row.input_rmse),
                row.species.clone(),
                format_float_trim(row.sps),
                row.true_phenotype.clone().unwrap_or_default(),
            ])?;
        } else {
            wtr.write_record([
                row.species_combo.clone(),
                format_float_trim(row.lambda1),
                format_float_trim(row.lambda2),
                format_float_trim(row.penalty_term),
                row.num_genes.to_string(),
                format_float_trim(row.input_rmse),
                row.species.clone(),
                format_float_trim(row.sps),
            ])?;
        }
    }

    wtr.flush()?;
    Ok(path)
}

fn write_gene_ranks_csv_single(
    output_dir: &Path,
    base_name: &str,
    genes: &[GeneAggregate],
    show_selected_sites: bool,
) -> Result<()> {
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by(|a, b| {
        let ar = genes[*a].single_best_rank.unwrap_or(usize::MAX);
        let br = genes[*b].single_best_rank.unwrap_or(usize::MAX);
        ar.cmp(&br).then_with(|| {
            genes[*b]
                .single_highest_gss
                .partial_cmp(&genes[*a].single_highest_gss)
                .unwrap_or(Ordering::Equal)
        })
    });

    let path = output_dir.join(format!("{}_gene_ranks.csv", base_name));
    let mut wtr =
        Writer::from_path(&path).with_context(|| format!("failed to create {}", path.display()))?;

    if show_selected_sites {
        wtr.write_record([
            "gene_name",
            "highest_gss",
            "best_rank",
            "num_selected_sites",
        ])?;
    } else {
        wtr.write_record(["gene_name", "highest_gss", "best_rank"])?;
    }

    for idx in order {
        let g = &genes[idx];
        if show_selected_sites {
            wtr.write_record([
                g.name.clone(),
                format_float_trim(g.single_highest_gss),
                g.single_best_rank
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "None".to_string()),
                g.selected_sites.len().to_string(),
            ])?;
        } else {
            wtr.write_record([
                g.name.clone(),
                format_float_trim(g.single_highest_gss),
                g.single_best_rank
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            ])?;
        }
    }

    wtr.flush()?;
    Ok(())
}

fn write_gene_ranks_csv_multimatrix(
    output_dir: &Path,
    base_name: &str,
    genes: &[GeneAggregate],
    show_selected_sites: bool,
) -> Result<()> {
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by(|a, b| {
        genes[*b]
            .num_combos_ranked
            .cmp(&genes[*a].num_combos_ranked)
            .then_with(|| {
                genes[*b]
                    .num_combos_ranked_top
                    .cmp(&genes[*a].num_combos_ranked_top)
            })
            .then_with(|| {
                let ar = genes[*a].best_ever_rank.unwrap_or(usize::MAX);
                let br = genes[*b].best_ever_rank.unwrap_or(usize::MAX);
                ar.cmp(&br)
            })
            .then_with(|| {
                genes[*b]
                    .highest_ever_gss
                    .partial_cmp(&genes[*a].highest_ever_gss)
                    .unwrap_or(Ordering::Equal)
            })
    });

    let path = output_dir.join(format!("{}_gene_ranks.csv", base_name));
    let mut wtr =
        Writer::from_path(&path).with_context(|| format!("failed to create {}", path.display()))?;

    if show_selected_sites {
        wtr.write_record([
            "gene_name",
            "num_combos_ranked",
            "num_combos_ranked_top",
            "highest_ever_gss",
            "best_ever_rank",
            "num_selected_sites",
        ])?;
    } else {
        wtr.write_record([
            "gene_name",
            "num_combos_ranked",
            "num_combos_ranked_top",
            "highest_ever_gss",
            "best_ever_rank",
        ])?;
    }

    for idx in order {
        let g = &genes[idx];
        if show_selected_sites {
            wtr.write_record([
                g.name.clone(),
                g.num_combos_ranked.to_string(),
                g.num_combos_ranked_top.to_string(),
                format_float_trim(g.highest_ever_gss),
                g.best_ever_rank
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "None".to_string()),
                g.selected_sites.len().to_string(),
            ])?;
        } else {
            wtr.write_record([
                g.name.clone(),
                g.num_combos_ranked.to_string(),
                g.num_combos_ranked_top.to_string(),
                format_float_trim(g.highest_ever_gss),
                g.best_ever_rank
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            ])?;
        }
    }

    wtr.flush()?;
    Ok(())
}

fn write_selected_sites_csv(
    output_dir: &Path,
    base_name: &str,
    genes: &[GeneAggregate],
) -> Result<()> {
    let path = output_dir.join(format!("{}_selected_sites.csv", base_name));
    let mut wtr =
        Writer::from_path(&path).with_context(|| format!("failed to create {}", path.display()))?;

    wtr.write_record(["gene_name", "position", "pss"])?;

    let mut rows: Vec<(String, usize, f64)> = Vec::new();
    for g in genes {
        for (pos, pss) in &g.selected_sites {
            rows.push((g.name.clone(), *pos + 1, *pss));
        }
    }

    rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for (gene, pos, pss) in rows {
        wtr.write_record([gene, pos.to_string(), format_float_trim(pss)])?;
    }

    wtr.flush()?;
    Ok(())
}

fn format_float_trim(x: f64) -> String {
    let s = format!("{:.15}", x);
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

fn is_fasta_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|s| s.to_str()).map(|s| s.to_ascii_lowercase()),
        Some(ext) if ext == "fas" || ext == "fasta" || ext == "fa" || ext == "faa"
    )
}

fn build_lambda_grid(args: &Args) -> Result<Vec<(f64, f64)>> {
    if args.use_logspace {
        lambda_logspace_grid(
            args.initial_lambda1,
            args.final_lambda1,
            args.initial_lambda2,
            args.final_lambda2,
            args.num_log_points,
        )
    } else {
        lambda_linear_grid(
            args.initial_lambda1,
            args.final_lambda1,
            args.initial_lambda2,
            args.final_lambda2,
            args.lambda_step,
        )
    }
}

fn lambda_row_ranges(lambda_grid: &[(f64, f64)]) -> Vec<std::ops::Range<usize>> {
    let mut rows = Vec::new();
    if lambda_grid.is_empty() {
        return rows;
    }
    let mut start = 0usize;
    while start < lambda_grid.len() {
        let l1 = lambda_grid[start].0;
        let mut end = start + 1;
        while end < lambda_grid.len() && (lambda_grid[end].0 - l1).abs() <= 1e-12 {
            end += 1;
        }
        rows.push(start..end);
        start = end;
    }
    rows
}

fn solve_lambda_row(
    prep: &PreprocessedData,
    group_weights: &[f64],
    lambda_row: &[(f64, f64)],
    penalty: f64,
    continuous: bool,
    maxiter: usize,
    disable_ec: bool,
    lipschitz: f64,
    progress_counter: Option<Arc<AtomicUsize>>,
) -> Result<Vec<ModelResult>> {
    let mut out = Vec::with_capacity(lambda_row.len());

    for (lambda1, lambda2) in lambda_row {
        let mut result = if prep.features.is_empty() {
            intercept_only_model(prep, *lambda1, *lambda2, penalty, continuous)
        } else {
            let mut r = solve_sparse_group_lasso(
                &prep.x,
                &prep.y_model,
                &prep.features,
                &prep.genes,
                group_weights,
                *lambda1,
                *lambda2,
                continuous,
                maxiter,
                disable_ec,
                1e-4,
                lipschitz,
                &[],
                0.0,
            )?;
            r.penalty_term = penalty;
            r
        };
        if prep.features.is_empty() {
            result.penalty_term = penalty;
        }
        out.push(result);
        if let Some(counter) = &progress_counter {
            counter.fetch_add(1, AtomicOrdering::Relaxed);
        }
    }

    Ok(out)
}

fn build_penalty_terms(args: &Args, genes: &[GeneMeta]) -> Result<Vec<f64>> {
    let kind = args.group_penalty_type.to_lowercase();
    if kind == "median" {
        let mut vars: Vec<usize> = genes
            .iter()
            .map(|g| g.var_site_count)
            .filter(|v| *v > 0)
            .collect();
        if vars.is_empty() {
            vars.push(1);
        }
        vars.sort_unstable();
        let median = if vars.len() % 2 == 1 {
            vars[vars.len() / 2] as f64
        } else {
            let hi = vars.len() / 2;
            let lo = hi - 1;
            ((vars[lo] + vars[hi]) as f64) / 2.0
        };
        return Ok(vec![median.floor()]);
    }

    linear_space(args.initial_gp_value, args.final_gp_value, args.gp_step)
}

fn median_var_sites_from_alignments(alignments: &[GeneAlignment]) -> Result<f64> {
    let mut vars: Vec<usize> = Vec::new();
    for gene in alignments {
        let seqs: Vec<Vec<u8>> = gene.seqs.values().cloned().collect();
        if seqs.is_empty() {
            continue;
        }
        let num_var = count_var_sites_python(&seqs);
        if num_var > 0 {
            vars.push(num_var);
        }
    }
    if vars.is_empty() {
        bail!("no alignments with variable sites found for median group penalty calculation");
    }
    vars.sort_unstable();
    let median = if vars.len() % 2 == 1 {
        vars[vars.len() / 2] as f64
    } else {
        let hi = vars.len() / 2;
        let lo = hi - 1;
        ((vars[lo] + vars[hi]) as f64) / 2.0
    };
    Ok(median.floor())
}

fn linear_space(start: f64, end: f64, step: f64) -> Result<Vec<f64>> {
    if step <= 0.0 {
        bail!("step must be > 0");
    }
    if end + 1e-12 < start {
        bail!("range end must be >= start");
    }

    let mut out = Vec::new();
    let mut val = start;
    while val <= end + 1e-12 {
        out.push(round_to(val, 6));
        val += step;
    }
    Ok(out)
}

fn logspace(start: f64, end: f64, points: usize) -> Result<Vec<f64>> {
    if start <= 0.0 || end <= 0.0 {
        bail!("logspace bounds must be > 0");
    }
    if points < 2 {
        bail!("logspace requires at least 2 points");
    }

    let log_start = start.log10();
    let log_end = end.log10();
    let step = (log_end - log_start) / (points as f64 - 1.0);

    let mut out = Vec::with_capacity(points);
    for i in 0..points {
        out.push(10f64.powf(log_start + (i as f64) * step));
    }
    Ok(out)
}

fn lambda_linear_grid(
    initial_lambda1: f64,
    final_lambda1: f64,
    initial_lambda2: f64,
    final_lambda2: f64,
    lambda_step: f64,
) -> Result<Vec<(f64, f64)>> {
    if lambda_step <= 0.0 {
        bail!("lambda_step must be > 0");
    }
    if final_lambda1 + 1e-12 < initial_lambda1 || final_lambda2 + 1e-12 < initial_lambda2 {
        bail!("final lambda must be >= initial lambda");
    }

    let mut out = Vec::new();
    let mut l1 = initial_lambda1;
    while l1 <= final_lambda1 + 1e-12 {
        let mut l2 = initial_lambda2;
        while l2 <= final_lambda2 + 1e-12 {
            out.push((l1, l2));
            l2 = round_to(l2 + lambda_step, 3);
        }
        l1 = round_to(l1 + lambda_step, 3);
    }
    Ok(out)
}

fn lambda_logspace_grid(
    initial_lambda1: f64,
    final_lambda1: f64,
    initial_lambda2: f64,
    final_lambda2: f64,
    num_values: usize,
) -> Result<Vec<(f64, f64)>> {
    let lambda1_values = logspace(initial_lambda1, final_lambda1, num_values)?;
    let lambda2_values = logspace(initial_lambda2, final_lambda2, num_values)?;
    let digits_to_round = initial_lambda1.log10().abs() as i32 + 5;
    let mut out = Vec::with_capacity(lambda1_values.len() * lambda2_values.len());
    for l1 in lambda1_values {
        for l2 in &lambda2_values {
            out.push((
                round_to(l1, digits_to_round as usize),
                round_to(*l2, digits_to_round as usize),
            ));
        }
    }
    Ok(out)
}

fn round_to(x: f64, digits: usize) -> f64 {
    let p = 10f64.powi(digits as i32);
    (x * p).round() / p
}

fn compute_group_weights(kind: &str, penalty_term: f64, genes: &[GeneMeta]) -> Vec<f64> {
    let k = kind.to_lowercase();
    genes
        .iter()
        .map(|g| {
            let feature_len = if g.feature_end >= g.feature_start {
                g.feature_end - g.feature_start + 1
            } else {
                0
            };
            match k.as_str() {
                "std" => (feature_len as f64).sqrt(),
                "linear" | "median" => g.var_site_count as f64 + penalty_term,
                "sqrt" => (g.var_site_count as f64).sqrt(),
                _ => (feature_len as f64).sqrt(),
            }
        })
        .collect()
}

fn initial_intercept(y: &[f64], continuous: bool) -> f64 {
    if y.is_empty() {
        return 0.0;
    }
    if continuous {
        return y.iter().sum::<f64>() / y.len() as f64;
    }
    let pos = y.iter().filter(|v| **v > 0.0).count() as f64;
    let neg = y.len() as f64 - pos;
    if pos <= 0.0 || neg <= 0.0 {
        return 0.0;
    }
    (pos / neg).ln()
}

#[allow(clippy::too_many_arguments)]
fn solve_sparse_group_lasso(
    x: &[Vec<f64>],
    y: &[f64],
    features: &[FeatureMeta],
    genes: &[GeneMeta],
    group_weights: &[f64],
    lambda1_frac: f64,
    lambda2_frac: f64,
    continuous: bool,
    maxiter: usize,
    disable_ec: bool,
    _tol: f64,
    _lipschitz: f64,
    _warm_beta: &[f64],
    _warm_intercept: f64,
) -> Result<ModelResult> {
    let n = x.len();
    let p = features.len();
    if n == 0 || p == 0 {
        bail!("invalid matrix shape for solver");
    }

    let mut beta = vec![0.0_f64; p];
    let intercept: f64;
    let scores: Vec<f64>;
    let r_flag_scaled = true;

    if continuous {
        let aty = mat_t_vec_mul(x, y);
        let lambda1_abs;
        let lambda2_abs;
        if r_flag_scaled {
            let lambda1_max = aty.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            lambda1_abs = lambda1_frac * lambda1_max;
            let mut tmp = vec![0.0_f64; p];
            for j in 0..p {
                tmp[j] = (aty[j].abs() - lambda1_abs).max(0.0);
            }
            let lambda2_max = compute_lambda2_max(&tmp, genes, group_weights);
            lambda2_abs = lambda2_frac * lambda2_max;
        } else {
            lambda1_abs = lambda1_frac;
            lambda2_abs = lambda2_frac;
        }

        let mut ax = vec![0.0_f64; n];
        let mut xp = beta.clone();
        let mut axp = ax.clone();
        let mut xxp = vec![0.0_f64; p];
        let mut alphap = 0.0_f64;
        let mut alpha = 1.0_f64;
        let mut l_const = 1.0_f64;
        let mut bflag = false;

        let mut s = vec![0.0_f64; p];
        let mut as_vec = vec![0.0_f64; n];
        let mut atas = vec![0.0_f64; p];
        let mut g = vec![0.0_f64; p];
        let mut v = vec![0.0_f64; p];
        let mut dv = vec![0.0_f64; p];
        let mut av = vec![0.0_f64; n];

        for iter_step in 0..maxiter {
            let beta_fac = (alphap - 1.0) / alpha;
            for j in 0..p {
                s[j] = beta[j] + xxp[j] * beta_fac;
            }
            for i in 0..n {
                as_vec[i] = ax[i] + (ax[i] - axp[i]) * beta_fac;
            }

            mat_t_vec_mul_into(x, &as_vec, &mut atas);
            for j in 0..p {
                g[j] = atas[j] - aty[j];
            }
            xp.clone_from(&beta);
            axp.clone_from(&ax);

            loop {
                for j in 0..p {
                    v[j] = s[j] - g[j] / l_const;
                }
                beta = altra_vector(
                    &v,
                    lambda1_abs / l_const,
                    lambda2_abs / l_const,
                    genes,
                    group_weights,
                );

                for j in 0..p {
                    dv[j] = beta[j] - s[j];
                }
                mat_vec_mul_into(x, &beta, &mut ax);

                for i in 0..n {
                    av[i] = ax[i] - as_vec[i];
                }
                let r_sum = dot(&dv, &dv);
                let l_sum = dot(&av, &av);

                if r_sum <= 1e-20 {
                    bflag = true;
                    break;
                }

                let target = r_sum * l_const;
                if line_search_accept_legacy_f64(l_sum, target, disable_ec) {
                    break;
                }
                l_const = (2.0 * l_const).max(l_sum / r_sum);
            }

            alphap = alpha;
            alpha = (1.0 + (4.0 * alpha * alpha + 1.0).sqrt()) / 2.0;

            for j in 0..p {
                xxp[j] = beta[j] - xp[j];
            }
            let _tree_norm = tree_norm(&beta, lambda1_abs, lambda2_abs, genes, group_weights);
            if bflag {
                break;
            }
            if (iter_step + 1) % maxiter == 0 {
                alphap = 0.0;
                alpha = 1.0;
                xp.clone_from(&beta);
                axp.clone_from(&ax);
                xxp.fill(0.0);
                l_const /= 2.0;
            }
        }

        scores = ax;
        intercept = 0.0;
    } else {
        let sample_weights = vec![1.0_f64 / n as f64; n];
        let p_flag: Vec<usize> = y
            .iter()
            .enumerate()
            .filter_map(|(i, v)| {
                if (*v - 1.0).abs() < 1e-12 {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
        let not_p_flag: Vec<usize> = y
            .iter()
            .enumerate()
            .filter_map(|(i, v)| {
                if (*v - 1.0).abs() < 1e-12 {
                    None
                } else {
                    Some(i)
                }
            })
            .collect();

        let sum_w: f64 = sample_weights.iter().sum();
        let m1 = p_flag.iter().map(|i| sample_weights[*i]).sum::<f64>() / sum_w;
        let m2 = 1.0_f64 - m1;

        let mut b = vec![0.0_f64; n];
        for i in &p_flag {
            b[*i] = m2;
        }
        for i in &not_p_flag {
            b[*i] = -m1;
        }
        for i in 0..n {
            b[i] *= sample_weights[i];
        }

        let atb = mat_t_vec_mul(x, &b);
        let lambda1_max = atb.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let lambda1_abs = lambda1_frac * lambda1_max;
        let mut tmp = vec![0.0_f64; p];
        for j in 0..p {
            tmp[j] = (atb[j].abs() - lambda1_abs).max(0.0);
        }
        let lambda2_max = compute_lambda2_max(&tmp, genes, group_weights);
        let lambda2_abs = lambda2_frac * lambda2_max;

        let mut intercept64 = if m2 > 0.0 { (m1 / m2).ln() } else { 0.0 };
        let mut beta64 = vec![0.0_f64; p];
        let mut ax = vec![0.0_f64; n];
        let mut xp = beta64.clone();
        let mut axp = ax.clone();
        let mut xxp = vec![0.0_f64; p];
        let mut ccp = 0.0_f64;
        let mut alphap = 0.0_f64;
        let mut alpha = 1.0_f64;
        let mut l_const = 1.0_f64 / n as f64;
        let mut bflag = false;

        let mut weighty = vec![0.0_f64; n];
        for i in 0..n {
            weighty[i] = sample_weights[i] * y[i];
        }

        let mut s = vec![0.0_f64; p];
        let mut as_vec = vec![0.0_f64; n];
        let mut aa = vec![0.0_f64; n];
        let mut bb = vec![0.0_f64; n];
        let mut prob = vec![0.0_f64; n];
        let mut g = vec![0.0_f64; p];
        let mut v = vec![0.0_f64; p];
        let mut dv = vec![0.0_f64; p];

        for iter_step in 0..maxiter {
            let beta_fac = (alphap - 1.0_f64) / alpha;
            for j in 0..p {
                s[j] = beta64[j] + xxp[j] * beta_fac;
            }
            let sc = intercept64 + beta_fac * ccp;

            for i in 0..n {
                as_vec[i] = ax[i] + (ax[i] - axp[i]) * beta_fac;
            }

            for i in 0..n {
                aa[i] = -y[i] * (as_vec[i] + sc);
                bb[i] = aa[i].max(0.0_f64);
            }
            let mut fun_s = 0.0_f64;
            for i in 0..n {
                fun_s +=
                    sample_weights[i] * (((-bb[i]).exp() + (aa[i] - bb[i]).exp()).ln() + bb[i]);
            }

            for i in 0..n {
                prob[i] = 1.0_f64 / (1.0_f64 + aa[i].exp());
                b[i] = -weighty[i] * (1.0_f64 - prob[i]);
            }
            let gc: f64 = b.iter().sum();
            mat_t_vec_mul_into(x, &b, &mut g);

            xp.clone_from(&beta64);
            axp.clone_from(&ax);
            let cp = intercept64;

            loop {
                for j in 0..p {
                    v[j] = s[j] - g[j] / l_const;
                }
                intercept64 = sc - gc / l_const;
                beta64 = altra_vector(
                    &v,
                    lambda1_abs / l_const,
                    lambda2_abs / l_const,
                    genes,
                    group_weights,
                );

                for j in 0..p {
                    dv[j] = beta64[j] - s[j];
                }

                mat_vec_mul_into(x, &beta64, &mut ax);
                for i in 0..n {
                    aa[i] = -y[i] * (ax[i] + intercept64);
                    bb[i] = aa[i].max(0.0_f64);
                }
                let mut fun_x = 0.0_f64;
                for i in 0..n {
                    fun_x +=
                        sample_weights[i] * (((-bb[i]).exp() + (aa[i] - bb[i]).exp()).ln() + bb[i]);
                }

                let r_sum = (dot(&dv, &dv) + (intercept64 - sc) * (intercept64 - sc)) / 2.0_f64;
                let l_sum = fun_x - fun_s - dot(&dv, &g) - (intercept64 - sc) * gc;

                if r_sum <= 1e-20_f64 {
                    bflag = true;
                    break;
                }
                let target = r_sum * l_const;
                if line_search_accept_legacy_f64(l_sum, target, disable_ec) {
                    break;
                }
                l_const = (2.0_f64 * l_const).max(l_sum / r_sum);
            }

            alphap = alpha;
            alpha = (1.0_f64 + (4.0_f64 * alpha * alpha + 1.0_f64).sqrt()) / 2.0_f64;

            for j in 0..p {
                xxp[j] = beta64[j] - xp[j];
            }
            ccp = intercept64 - cp;
            let _tree_norm = tree_norm(&beta64, lambda1_abs, lambda2_abs, genes, group_weights);
            if bflag {
                break;
            }
            if (iter_step + 1) % maxiter == 0 {
                alphap = 0.0_f64;
                alpha = 1.0_f64;
                xp.clone_from(&beta64);
                axp.clone_from(&ax);
                xxp.fill(0.0_f64);
                l_const /= 2.0_f64;
            }
        }

        beta = beta64;
        intercept = intercept64;
        let mut s = vec![0.0_f64; n];
        for i in 0..n {
            s[i] = ax[i] + intercept64;
        }
        scores = s;
    }

    let rmse = compute_rmse(&scores, y);

    let mut gene_gss = vec![0.0_f64; genes.len()];
    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let mut s = 0.0_f64;
        for b in &beta[g.feature_start..=g.feature_end] {
            s += b.abs();
        }
        gene_gss[gidx] = s;
    }

    let num_genes = gene_gss.iter().filter(|x| **x > 0.0).count();
    let selected_sites = collect_selected_sites(features, &beta);

    Ok(ModelResult {
        lambda1: lambda1_frac,
        lambda2: lambda2_frac,
        penalty_term: 0.0,
        beta,
        intercept,
        rmse,
        num_genes,
        gene_gss,
        selected_sites,
    })
}

fn collect_selected_sites(features: &[FeatureMeta], beta: &[f64]) -> Vec<SelectedSite> {
    let mut sums: HashMap<(usize, usize), f64> = HashMap::new();
    for (fm, w) in features.iter().zip(beta.iter()) {
        let p = w.abs();
        if p == 0.0 {
            continue;
        }
        let key = (fm.gene_idx, fm.position);
        sums.entry(key).and_modify(|v| *v += p).or_insert(p);
    }

    let mut out: Vec<SelectedSite> = sums
        .into_iter()
        .map(|((gene_idx, position), pss)| SelectedSite {
            gene_idx,
            position,
            pss,
        })
        .collect();
    out.sort_by(|a, b| {
        a.gene_idx
            .cmp(&b.gene_idx)
            .then(a.position.cmp(&b.position))
    });
    out
}

fn line_search_accept_legacy_f64(l_sum: f64, target: f64, disable_ec: bool) -> bool {
    if disable_ec {
        return l_sum <= target;
    }
    // Compatibility with legacy C++ behavior:
    // `abs(diff) < 1e-12` was evaluated with integer-style `abs`, effectively
    // treating any diff with truncated int value 0 as acceptable.
    let diff = l_sum - target;
    l_sum < target || (diff as i32).abs() == 0
}

fn line_search_accept_legacy_f32(l_sum: f32, target: f32, disable_ec: bool) -> bool {
    if disable_ec {
        return l_sum <= target;
    }
    let diff = l_sum - target;
    l_sum < target || (diff as i32).abs() == 0
}

#[cfg(test)]
mod tests {
    use super::{line_search_accept_legacy_f32, line_search_accept_legacy_f64};

    #[test]
    fn epsilon_comparison_accepts_small_positive_diff() {
        assert!(line_search_accept_legacy_f64(1.2, 1.0, false));
        assert!(line_search_accept_legacy_f32(1.2, 1.0, false));
    }

    #[test]
    fn disable_ec_restores_strict_acceptance() {
        assert!(!line_search_accept_legacy_f64(1.2, 1.0, true));
        assert!(!line_search_accept_legacy_f32(1.2, 1.0, true));
        assert!(line_search_accept_legacy_f64(1.0, 1.0, true));
        assert!(line_search_accept_legacy_f32(1.0, 1.0, true));
    }
}

fn predict_scores(x: &[Vec<f64>], beta: &[f64], intercept: f64) -> Vec<f64> {
    let mut out = vec![0.0_f64; x.len()];
    for (i, row) in x.iter().enumerate() {
        let mut s = intercept;
        for (v, b) in row.iter().zip(beta.iter()) {
            s += v * b;
        }
        out[i] = s;
    }
    out
}

fn compute_rmse(pred: &[f64], y: &[f64]) -> f64 {
    let n = pred.len().max(1) as f64;
    let mut s = 0.0_f64;
    for (a, b) in pred.iter().zip(y.iter()) {
        let d = a - b;
        s += d * d;
    }
    (s / n).sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn mat_vec_mul_into(x: &[Vec<f64>], beta: &[f64], out: &mut [f64]) {
    debug_assert_eq!(x.len(), out.len());
    for (i, row) in x.iter().enumerate() {
        let mut s = 0.0_f64;
        for (v, b) in row.iter().zip(beta.iter()) {
            s += v * b;
        }
        out[i] = s;
    }
}

fn mat_t_vec_mul_into(x: &[Vec<f64>], v: &[f64], out: &mut [f64]) {
    if x.is_empty() {
        return;
    }
    debug_assert_eq!(x.len(), v.len());
    debug_assert_eq!(x[0].len(), out.len());
    out.fill(0.0_f64);
    for (i, row) in x.iter().enumerate() {
        let vi = v[i];
        for (j, xij) in row.iter().enumerate() {
            out[j] += xij * vi;
        }
    }
}

fn mat_t_vec_mul(x: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    if x.is_empty() {
        return Vec::new();
    }
    let p = x[0].len();
    let mut out = vec![0.0_f64; p];
    mat_t_vec_mul_into(x, v, &mut out);
    out
}

fn dot_f32(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn mat_vec_mul_f32(x: &[Vec<f64>], beta: &[f32]) -> Vec<f32> {
    let mut out = vec![0.0_f32; x.len()];
    for (i, row) in x.iter().enumerate() {
        let mut s = 0.0_f32;
        for (v, b) in row.iter().zip(beta.iter()) {
            s += (*v as f32) * *b;
        }
        out[i] = s;
    }
    out
}

fn mat_t_vec_mul_f32(x: &[Vec<f64>], v: &[f32]) -> Vec<f32> {
    if x.is_empty() {
        return Vec::new();
    }
    let p = x[0].len();
    let mut out = vec![0.0_f32; p];
    for (i, row) in x.iter().enumerate() {
        let vi = v[i];
        for (j, xij) in row.iter().enumerate() {
            out[j] += (*xij as f32) * vi;
        }
    }
    out
}

fn altra_vector(
    v: &[f64],
    lambda1: f64,
    lambda2: f64,
    genes: &[GeneMeta],
    group_weights: &[f64],
) -> Vec<f64> {
    let mut x = vec![0.0_f64; v.len()];
    for j in 0..v.len() {
        if v[j] > lambda1 {
            x[j] = v[j] - lambda1;
        } else if v[j] < -lambda1 {
            x[j] = v[j] + lambda1;
        } else {
            x[j] = 0.0;
        }
    }

    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f64;
        for j in start..=end {
            two_norm += x[j] * x[j];
        }
        two_norm = two_norm.sqrt();
        let lambda = lambda2 * group_weights[gidx];
        if two_norm > lambda {
            let ratio = (two_norm - lambda) / two_norm;
            for j in start..=end {
                x[j] *= ratio;
            }
        } else {
            for j in start..=end {
                x[j] = 0.0;
            }
        }
    }
    x
}

fn tree_norm(
    x: &[f64],
    lambda1: f64,
    lambda2: f64,
    genes: &[GeneMeta],
    group_weights: &[f64],
) -> f64 {
    let mut t = lambda1 * x.iter().map(|v| v.abs()).sum::<f64>();
    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f64;
        for j in start..=end {
            two_norm += x[j] * x[j];
        }
        t += lambda2 * group_weights[gidx] * two_norm.sqrt();
    }
    t
}

fn altra_vector_f32(
    v: &[f32],
    lambda1: f32,
    lambda2: f32,
    genes: &[GeneMeta],
    group_weights: &[f32],
) -> Vec<f32> {
    let mut x = vec![0.0_f32; v.len()];
    for j in 0..v.len() {
        if v[j] > lambda1 {
            x[j] = v[j] - lambda1;
        } else if v[j] < -lambda1 {
            x[j] = v[j] + lambda1;
        } else {
            x[j] = 0.0;
        }
    }

    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f32;
        for j in start..=end {
            two_norm += x[j] * x[j];
        }
        two_norm = two_norm.sqrt();
        let lambda = lambda2 * group_weights[gidx];
        if two_norm > lambda {
            let ratio = (two_norm - lambda) / two_norm;
            for j in start..=end {
                x[j] *= ratio;
            }
        } else {
            for j in start..=end {
                x[j] = 0.0;
            }
        }
    }
    x
}

fn tree_norm_f32(
    x: &[f32],
    lambda1: f32,
    lambda2: f32,
    genes: &[GeneMeta],
    group_weights: &[f32],
) -> f32 {
    let mut t = lambda1 * x.iter().map(|v| v.abs()).sum::<f32>();
    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f32;
        for j in start..=end {
            two_norm += x[j] * x[j];
        }
        t += lambda2 * group_weights[gidx] * two_norm.sqrt();
    }
    t
}

fn compute_lambda2_max(temp: &[f64], genes: &[GeneMeta], group_weights: &[f64]) -> f64 {
    let mut lambda2_max = 0.0_f64;
    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let denom = group_weights[gidx];
        if denom <= 0.0 {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f64;
        for j in start..=end {
            two_norm += temp[j] * temp[j];
        }
        let val = two_norm.sqrt() / denom;
        if val > lambda2_max {
            lambda2_max = val;
        }
    }
    lambda2_max
}

fn compute_lambda2_max_f32(temp: &[f32], genes: &[GeneMeta], group_weights: &[f32]) -> f32 {
    let mut lambda2_max = 0.0_f32;
    for (gidx, g) in genes.iter().enumerate() {
        if g.feature_end < g.feature_start {
            continue;
        }
        let denom = group_weights[gidx];
        if denom <= 0.0 {
            continue;
        }
        let start = g.feature_start;
        let end = g.feature_end;
        let mut two_norm = 0.0_f32;
        for j in start..=end {
            two_norm += temp[j] * temp[j];
        }
        let val = two_norm.sqrt() / denom;
        if val > lambda2_max {
            lambda2_max = val;
        }
    }
    lambda2_max
}

fn estimate_lipschitz(x: &[Vec<f64>], continuous: bool) -> f64 {
    let n = x.len();
    if n == 0 {
        return 1.0;
    }
    let p = x[0].len();
    if p == 0 {
        return 1.0;
    }

    let mut v = vec![1.0_f64 / (p as f64).sqrt(); p];
    let mut lambda = 1.0_f64;

    for _ in 0..20 {
        let mut u = vec![0.0_f64; n];
        for i in 0..n {
            let mut s = 0.0_f64;
            for (xij, vj) in x[i].iter().zip(v.iter()) {
                s += xij * vj;
            }
            u[i] = s;
        }

        let mut v2 = vec![0.0_f64; p];
        for i in 0..n {
            let ui = u[i];
            for (j, xij) in x[i].iter().enumerate() {
                v2[j] += xij * ui;
            }
        }

        let norm = v2.iter().map(|z| z * z).sum::<f64>().sqrt();
        if norm <= 1e-12 {
            lambda = 1.0;
            break;
        }

        for j in 0..p {
            v[j] = v2[j] / norm;
        }

        lambda = v.iter().zip(v2.iter()).map(|(a, b)| a * b).sum::<f64>();
    }

    let mut l = lambda / (n as f64);
    if !continuous {
        l *= 0.25;
    }
    l.max(1e-6)
}

fn python_round5_tag(x: f64) -> String {
    let rounded = (x * 100000.0).round() / 100000.0;
    let mut s = rounded.to_string();
    if !s.contains('.') && !s.contains('e') && !s.contains('E') {
        s.push_str(".0");
    }
    s
}

fn write_model_file(
    model_dir: &Path,
    base: &str,
    features: &[FeatureMeta],
    run: &ModelResult,
    combo_label: Option<&str>,
) -> Result<()> {
    let l1 = python_round5_tag(run.lambda1);
    let l2 = python_round5_tag(run.lambda2);

    let fname = if let Some(combo) = combo_label {
        format!(
            "{}_{}_l1_{}_l2_{}_out_feature_weights.txt",
            base, combo, l1, l2
        )
    } else {
        format!("{}_l1_{}_l2_{}_out_feature_weights.txt", base, l1, l2)
    };

    let path = model_dir.join(fname);
    let mut f = File::create(&path)
        .with_context(|| format!("failed to create model file {}", path.display()))?;

    for (fm, w) in features.iter().zip(run.beta.iter()) {
        if *w == 0.0 {
            continue;
        }
        writeln!(f, "{}\t{:.17e}", fm.label, w)?;
    }
    writeln!(f, "Intercept\t{:.17e}", run.intercept)?;
    Ok(())
}
