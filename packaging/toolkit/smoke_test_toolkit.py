#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
import zipfile
from pathlib import Path


def extract_archive(artifact: Path, out_dir: Path) -> None:
    suffixes = "".join(artifact.suffixes[-2:])
    if suffixes == ".tar.gz":
        with tarfile.open(artifact, "r:gz") as handle:
            handle.extractall(out_dir)
        return
    if artifact.suffix.lower() == ".zip":
        with zipfile.ZipFile(artifact, "r") as handle:
            handle.extractall(out_dir)
        return
    raise ValueError(f"unsupported artifact type: {artifact.name}")


def find_toolkit_root(extract_dir: Path) -> Path:
    candidates = [extract_dir / "esl-psc-toolkit"]
    candidates.extend(p for p in extract_dir.iterdir() if p.is_dir())
    for candidate in candidates:
        if (candidate / "bin").is_dir() and (candidate / "requirements-toolkit.txt").is_file():
            return candidate
    raise FileNotFoundError("unable to locate toolkit root in extracted artifact")


def run_checked(command: list[str], timeout: int = 60) -> None:
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}"
        )


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum_and_manifest(artifact: Path) -> None:
    checksum_file = artifact.parent / f"{artifact.name}.sha256"
    manifest_file = artifact.parent / f"{artifact.name}.manifest.json"
    if not checksum_file.is_file():
        raise FileNotFoundError(f"checksum file not found: {checksum_file}")
    if not manifest_file.is_file():
        raise FileNotFoundError(f"manifest file not found: {manifest_file}")

    raw = checksum_file.read_text(encoding="utf-8").strip()
    pieces = raw.split()
    if len(pieces) < 2:
        raise ValueError(f"invalid checksum file format: {checksum_file}")
    expected_hash = pieces[0]
    expected_name = pieces[-1].lstrip("*")
    if expected_name != artifact.name:
        raise ValueError(
            f"checksum references {expected_name}, expected {artifact.name}"
        )

    actual_hash = sha256sum(artifact)
    if actual_hash != expected_hash:
        raise ValueError(f"checksum mismatch for {artifact.name}")

    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("archive_file") != artifact.name:
        raise ValueError("manifest archive_file mismatch")
    if manifest.get("archive_sha256") != actual_hash:
        raise ValueError("manifest archive_sha256 mismatch")
    if int(manifest.get("archive_size_bytes", -1)) != artifact.stat().st_size:
        raise ValueError("manifest archive_size_bytes mismatch")


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke test an ESL-PSC toolkit archive.")
    parser.add_argument("--artifact", required=True, help="Path to toolkit archive (.tar.gz or .zip)")
    args = parser.parse_args()

    artifact = Path(args.artifact).resolve()
    if not artifact.is_file():
        raise FileNotFoundError(f"artifact not found: {artifact}")
    verify_checksum_and_manifest(artifact)

    with tempfile.TemporaryDirectory(prefix="esl-psc-toolkit-smoke-") as temp_root:
        extract_dir = Path(temp_root) / "extract"
        extract_dir.mkdir(parents=True, exist_ok=True)
        extract_archive(artifact, extract_dir)
        toolkit_root = find_toolkit_root(extract_dir)
        bin_dir = toolkit_root / "bin"

        if artifact.suffix.lower() == ".zip":
            esl_main = str(bin_dir / "esl-psc.exe")
            run_checked([esl_main, "--help"])
            run_checked([esl_main, "pairs", "--help"])
            run_checked([esl_main, "site-counter", "--help"])
            run_checked([esl_main, "plot", "--help"])
        else:
            esl_main = str(bin_dir / "esl-psc")
            run_checked([esl_main, "--help"])
            run_checked([esl_main, "pairs", "--help"])
            run_checked([esl_main, "site-counter", "--help"])
            run_checked([esl_main, "plot", "--help"])

    print(f"smoke test passed: {artifact}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
