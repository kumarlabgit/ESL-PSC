#!/usr/bin/env bash
set -euo pipefail

cargo build --release --manifest-path esl_psc_rs/Cargo.toml

BIN_PATH=""
for candidate in \
  "esl_psc_rs/target/release/esl-psc" \
  "target/release/esl-psc"; do
  if [ -f "${candidate}" ]; then
    BIN_PATH="${candidate}"
    break
  fi
done
if [ -z "${BIN_PATH}" ]; then
  echo "could not find built binary esl-psc in expected target directories" >&2
  exit 1
fi

install -d "${PREFIX}/bin"
install -m 0755 "${BIN_PATH}" "${PREFIX}/bin/esl-psc"

install -d "${SP_DIR}/esl_psc_cli"
cp -a esl_psc_cli/. "${SP_DIR}/esl_psc_cli/"

install -d "${SP_DIR}/gui/core"
install -m 0644 gui/__init__.py "${SP_DIR}/gui/__init__.py"
install -m 0644 gui/core/fast_scan.py "${SP_DIR}/gui/core/fast_scan.py"
install -m 0644 gui/core/fasta_io.py "${SP_DIR}/gui/core/fasta_io.py"
install -m 0644 gui/core/ancestral_reconstruction.py "${SP_DIR}/gui/core/ancestral_reconstruction.py"
