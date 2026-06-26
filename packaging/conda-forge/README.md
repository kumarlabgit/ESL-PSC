# Conda-forge Packaging

This directory contains a conda-forge-oriented recipe and helper tooling.

## Render recipe for a tag

```bash
python3 packaging/conda-forge/render_meta_from_tag.py \
  --repo John-Allard/ESL-PSC \
  --tag v2.5.3 \
  --output dist/conda-forge/meta.yaml
```

This computes the source tarball SHA256 from GitHub and writes a concrete
`meta.yaml` suitable for staged-recipes/feedstock updates.

## Build locally (after rendering)

```bash
conda build packaging/conda-forge/recipe
```

## Notes

- Recipe installs:
  - `esl-psc` binary
  - `esl_psc_cli` module
  - minimal `gui.core` modules used by utility subcommands
- Runtime dependencies mirror toolkit requirements.
