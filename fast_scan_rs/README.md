# Fast Scan Rust CLI

The `fast_scan_rs` crate provides a compiled backend for ESL-PSC's fast scan
workflow. It mirrors the logic used by the Python implementation but uses
Rust's parallel iterators to dramatically reduce scan time on large alignment
collections.

The executable consumes a JSON description of the scan on **stdin** and writes
its results as a JSON array on **stdout**. Progress updates (if enabled) are
emitted to **stderr** so they will not interfere with downstream tools that
expect machine-readable output.

## Building the binary

```bash
cd fast_scan_rs
cargo build --release
```

The optimized binary will be created at
`fast_scan_rs/target/release/fast_scan_rs`.  The Python CLI (`esl_psc_cli
fast_scan_cli.py`) automatically prefers this binary when it is present and
executable; otherwise it falls back to the pure-Python implementation.

## JSON input schema

The executable expects a single JSON object that matches the structure below.
Optional fields may be omitted. The values correspond to the parameters that
ESL-PSC's Python fast scan uses internally.

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
  "min_out_ctrl_agreement": 1.0
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

## Direct usage example

```bash
cat <<'JSON' | fast_scan_rs/target/release/fast_scan_rs > results.json
{
  "alignment_dir": "photosynthesis_alignments",
  "combos": [
    {"conv": ["Maize"], "ctrl": ["Sorghum"]}
  ],
  "outgroup": "Oryza_sativa",
  "emit_progress": true
}
JSON
```

The resulting `results.json` file contains the summarized statistics for each
alignment processed by the binary.

## Integration with the Python CLI

`python -m esl_psc_cli.fast_scan_cli` wraps this binary. When a compatible
`fast_scan_rs` build is available, the CLI streams the required JSON over stdin
and performs the same post-processing, CSV export, and ranking logic as the GUI
fast scan workflow. Set the environment variable `FAST_SCAN_RS_DISABLE=1` to
force the Python fallback if needed.
