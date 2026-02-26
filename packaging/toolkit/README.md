# ESL-PSC Toolkit

This package contains command-line tools for ESL-PSC:

- `esl-psc` - unified main ESL-PSC CLI
- `esl-psc pairs` (or `esl-psc-pairs`) - auto pair-selection CLI
- `esl-psc site-counter` (or `site-counter`) - site counter CLI
- `esl-psc plot` (or `esl-psc-plot`) - plotting helper CLI

The toolkit is designed to run with your system Python. It does not ship a
second embedded Python runtime.

## Install

1. Install Python dependencies:

```bash
python3 -m pip install -r requirements-toolkit.txt
```

2. Add `bin/` to your `PATH` (or run commands with explicit `bin/...` paths).

If you relocate the `python/` folder away from the default toolkit layout,
set `ESL_PSC_PYTHONPATH` to that folder before running subcommands like
`esl-psc pairs`, `esl-psc site-counter`, or `esl-psc plot`.

## Basic Usage

```bash
esl-psc --help
esl-psc pairs --help
esl-psc site-counter --help
esl-psc plot --help
```

If these binaries are not on your `PATH`, run them from `bin/` directly.
