# Conda Packaging (WIP)

This directory contains a local conda-build recipe for packaging the unified
`esl-psc` CLI plus required Python modules for utility subcommands.

## Local Build

From repository root:

```bash
conda-build packaging/conda/recipe
```

To place built packages in a local folder:

```bash
conda-build packaging/conda/recipe --output-folder dist/conda
```

## Notes

- The recipe installs:
  - `esl-psc` binary
  - `esl_psc_cli` Python package
  - minimal `gui.core` modules used by CLI utilities
- Runtime Python dependencies mirror `requirements-toolkit.txt`.
- For conda-forge recipe generation tooling, see `packaging/conda-forge/`.
