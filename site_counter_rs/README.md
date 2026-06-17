# Site Counter Backend

The `site_counter_rs` crate provides the compiled Site Counter backend for
ESL-PSC. Site Counter analyses should usually be run through the unified
command-line interface:

```bash
esl-psc site-counter --help
```

See [`docs/commands/site-counter.md`](../docs/commands/site-counter.md) for
the command documentation.

The backend itself consumes a JSON description of the scan on **stdin** and
writes its results as a JSON array on **stdout**. Progress updates, if enabled,
are emitted to **stderr** so they do not interfere with tools that expect
machine-readable output.

This JSON interface is used by ESL-PSC and is also useful for development and
packaging tests.

## Building the binary

```bash
cd site_counter_rs
cargo build --release
```

The optimized binary will be created at
`site_counter_rs/target/release/site_counter_rs`.

## JSON input schema

The executable expects a single JSON object that matches the structure below.
Optional fields may be omitted. The values correspond to the backend request
created by the GUI and `esl-psc site-counter` wrappers.

```json
{
  "alignment_dir": "path/to/alignment_dir",
  "files": ["Gene1.fas", "Gene2.fas"],
  "combos": [
    {"conv": ["SpeciesA", "SpeciesB"], "ctrl": ["SpeciesC", "SpeciesD"]}
  ],
  "outgroup": "OutgroupSpecies",
  "cs_threshold": 4,
  "emit_progress": true,
  "min_out_ctrl_agreement": 1.0,
  "tree_file": "path/to/tree.nwk",
  "tree_json": { "name": null, "children": [...] },
  "analysis_species": ["SpeciesA", "SpeciesB", "SpeciesC", "SpeciesD"]
}
```

* `alignment_dir` (string, required): directory containing FASTA alignment
  files. When `files` is omitted every `*.fas` file in the directory is
  scanned.
* `files` (list of strings, optional): explicit list of alignment filenames to
  analyze.
* `combos` (list of objects, required): each object defines the convergent
  (`conv`) and control (`ctrl`) species identifiers for one contrast pairing.
  The binary automatically filters out missing species before scoring each
  alignment.
* `outgroup` (string, required): the species identifier used for outgroup
  comparison.
* `cs_threshold` (integer, optional, default `4`): minimum raw CCS score for a
  site to be tallied as `cs_sites_ge_4`.
* `emit_progress` (boolean, optional, default `false`): when true, prints a
  `PROGRESS <count>` line to stderr every 200 files (and at completion).
* `min_out_ctrl_agreement` (float in `[0,1]`, optional, default `1.0`): minimum
  fraction of control residues that must match the outgroup residue before a
  site is counted as convergent. Values outside the range are clamped.
* `tree_file` (string, optional): path to a Newick or NEXUS tree file for
  parsimony-based ancestral reconstruction. When provided with `analysis_species`,
  the MRCA of analysis species is reconstructed and used as the outgroup.
* `tree_json` (object, optional): pre-parsed tree structure (from Python/Biopython)
  with schema `{"name": str | null, "children": [...]}`. Preferred over `tree_file`
  when both are present.
* `analysis_species` (list of strings, optional): species in the convergence analysis.
  Required when using ancestral reconstruction to identify the MRCA node.
* `require_unambiguous_mrca` (boolean, optional, default `false`): when true, sites
  where ancestral reconstruction yields multiple equally parsimonious residues are
  excluded from CCS detection. By default (false), such ambiguous sites are evaluated
  using the set of possible ancestral states: a site counts as convergent if controls
  match any possible ancestral state while convergent species share a derived residue
  not in the ancestral set.

## Output format

The program prints a JSON array. Each element summarizes one alignment file:

```json
{
  "gene": "Gene1",
  "avg_true": 1.25,
  "avg_control": 0.1,
  "diff": 1.15,
  "variable_sites": 42,
  "cs_sites_ge_4": 5,
  "k_pairs": 3,
  "per_combo_true": [2.0, null],
  "per_combo_diff": [1.5, null]
}
```

* `gene`: alignment filename without the `.fas` extension.
* `avg_true` / `avg_control`: mean CCS counts across eligible convergent and
  control pairings.
* `diff`: difference `avg_true - avg_control`.
* `variable_sites`: number of polymorphic positions that survive filtering.
* `cs_sites_ge_4`: highest count of convergent sites meeting `cs_threshold`
  across all combos for the alignment.
* `k_pairs`: number of convergent/control pairs that produced that maximum
  `cs_sites_ge_4` value.
* `per_combo_true`: per-combo CCS counts (or `null` when a combo was
  ineligible).
* `per_combo_diff`: per-combo difference between convergent and control CCS
  counts when both sides were eligible.

## Direct Backend Invocation

Direct backend use requires a JSON request file. For ordinary Site Counter
analyses, use `esl-psc site-counter` instead.

Example `request.json`:

```json
{
  "alignment_dir": "demo_data/photosynthesis/alignments",
  "combos": [
    {"conv": ["Maize"], "ctrl": ["Sorghum"]}
  ],
  "outgroup": "Oryza_sativa",
  "emit_progress": true
}
```

Run the standalone backend:

```bash
site_counter_rs/target/release/site_counter_rs < request.json > results.json
```

The same backend can also be reached through the unified executable:

```bash
esl-psc site-counter-backend < request.json > results.json
```

The resulting `results.json` file contains the summarized statistics for the
processed alignments.

## Integration with ESL-PSC

The `esl-psc site-counter` subcommand and the GUI Site Counter workflow call
this backend automatically. They handle species-groups parsing, CSV export,
ranking, and progress reporting.
