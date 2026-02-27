# APT / Debian Packaging

This directory contains tooling to build a Debian package for ESL-PSC.

## Build locally

From repository root:

```bash
cargo build --release --manifest-path esl_psc_rs/Cargo.toml
python3 packaging/apt/build_deb.py \
  --dest-dir dist/apt \
  --rust-main esl_psc_rs/target/release/esl-psc
```

Output:

- `dist/apt/esl-psc_<version>_amd64.deb`

## Install test

```bash
sudo apt install ./dist/apt/esl-psc_<version>_amd64.deb
esl-psc --help
esl-psc pairs --help
esl-psc site-counter --help
esl-psc plot --help
```

