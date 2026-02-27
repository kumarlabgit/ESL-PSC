# Homebrew Packaging

This directory contains tooling to generate a Homebrew formula from published
release assets.

## Generate formula from a release tag

```bash
python3 packaging/homebrew/render_formula_from_release.py \
  --repo John-Allard/ESL-PSC \
  --tag v2.4.1 \
  --output dist/homebrew/esl-psc.rb
```

The script reads release assets:

- `esl-psc-toolkit-v<version>-linux-x86_64.tar.gz`
- `esl-psc-toolkit-v<version>-macos-x86_64.tar.gz`
- corresponding `.sha256` files

and renders a formula that:

- installs `esl-psc`
- configures `ESL_PSC_PYTHONPATH`
- points `ESL_PSC_PYTHON` to Homebrew Python

