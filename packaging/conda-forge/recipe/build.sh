#!/usr/bin/env bash
set -euo pipefail

cargo build --release --manifest-path esl_psc_rs/Cargo.toml

install -d "${PREFIX}/bin"
install -m 0755 esl_psc_rs/target/release/esl-psc "${PREFIX}/bin/esl-psc"

install -d "${SP_DIR}/esl_psc_cli"
cp -a esl_psc_cli/. "${SP_DIR}/esl_psc_cli/"

install -d "${SP_DIR}/gui/core"
install -m 0644 gui/__init__.py "${SP_DIR}/gui/__init__.py"
install -m 0644 gui/core/fast_scan.py "${SP_DIR}/gui/core/fast_scan.py"
install -m 0644 gui/core/fasta_io.py "${SP_DIR}/gui/core/fasta_io.py"
install -m 0644 gui/core/ancestral_reconstruction.py "${SP_DIR}/gui/core/ancestral_reconstruction.py"

