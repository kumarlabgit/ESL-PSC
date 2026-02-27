# Conda Packaging (WIP)

This directory contains a first-pass conda recipe for packaging the unified
`esl-psc` CLI plus required Python modules for utility subcommands.

## Local Build

From repository root:

```bash
conda build packaging/conda/recipe
```

To place built packages in a local folder:

```bash
conda build packaging/conda/recipe --output-folder dist/conda
```

## Notes

- The recipe installs:
  - `esl-psc` binary
  - `esl_psc_cli` Python package
  - minimal `gui.core` modules used by CLI utilities
- Runtime Python dependencies mirror `requirements-toolkit.txt`.
- This recipe is intended as the base for conda-forge/Anaconda distribution.

