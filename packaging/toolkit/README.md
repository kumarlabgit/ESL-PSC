# ESL-PSC Toolkit

This package contains command-line tools for ESL-PSC:

- `esl-psc` - unified main ESL-PSC CLI
- `esl-psc pairs` - auto pair-selection CLI
- `esl-psc site-counter` - site counter CLI

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
`esl-psc pairs` or `esl-psc site-counter`.

## Basic Usage

```bash
esl-psc --help
esl-psc pairs --help
esl-psc site-counter --help
```

If these binaries are not on your `PATH`, run them from `bin/` directly.

## Release Artifact Format

Toolkit release jobs publish three files per platform:

- `esl-psc-toolkit-v<version>-<os>-<arch>.<tar.gz|zip>`
- `esl-psc-toolkit-v<version>-<os>-<arch>.<tar.gz|zip>.sha256`
- `esl-psc-toolkit-v<version>-<os>-<arch>.<tar.gz|zip>.manifest.json`

The `.sha256` file can be verified with standard tooling:

```bash
sha256sum -c esl-psc-toolkit-v<version>-<os>-<arch>.<tar.gz|zip>.sha256
```
