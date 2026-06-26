# Homebrew Packaging

This directory contains tooling to generate a Homebrew formula from a GitHub
source tag.

## Generate formula from a tag

```bash
python3 packaging/homebrew/render_formula_from_release.py \
  --repo John-Allard/ESL-PSC \
  --tag v2.5.3 \
  --output dist/homebrew/esl-psc.rb
```

The script computes the SHA256 of
`https://github.com/<owner>/<repo>/archive/refs/tags/<tag>.tar.gz` and renders
a formula that:

- builds the Rust CLI from source with Cargo
- installs `esl-psc`
- configures `ESL_PSC_PYTHONPATH`
- points `ESL_PSC_PYTHON` to Homebrew Python
