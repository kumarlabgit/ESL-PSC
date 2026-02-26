# ESL-PSC Toolkit

This package contains command-line tools for ESL-PSC:

- `esl-psc` - unified main ESL-PSC CLI
- `esl-psc-pairs` - auto pair-selection CLI
- `site-counter` - site counter CLI
- `esl-psc-plot` - plotting helper CLI used by `esl-psc` plot flags
- `site_counter_rs` - Rust backend used by `site-counter`

The toolkit is designed to run with your system Python. It does not ship a
second embedded Python runtime.

## Install

1. Install Python dependencies:

```bash
python3 -m pip install -r requirements-toolkit.txt
```

2. Add `bin/` to your `PATH` (or run commands with explicit `bin/...` paths).

## Basic Usage

```bash
esl-psc --help
esl-psc-pairs --help
site-counter --help
esl-psc-plot --help
```

If these binaries are not on your `PATH`, run them from `bin/` directly.
