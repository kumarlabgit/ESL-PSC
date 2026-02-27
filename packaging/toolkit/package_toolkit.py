#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from stage_toolkit import stage_toolkit


def detect_version(repo_root: Path) -> str:
    cargo_toml = repo_root / "esl_psc_rs" / "Cargo.toml"
    in_package = False
    for raw_line in cargo_toml.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            in_package = line == "[package]"
            continue
        if in_package and line.startswith("version"):
            _, value = line.split("=", 1)
            return value.strip().strip('"').strip("'")
    raise ValueError(f"unable to detect version from {cargo_toml}")


def sanitize_version(version: str) -> str:
    cleaned = version.strip().replace("/", "-")
    if not cleaned:
        raise ValueError("version must be non-empty")
    return cleaned


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_tar_gz(source_dir: Path, archive_path: Path) -> None:
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(source_dir, arcname=source_dir.name)


def make_zip(source_dir: Path, archive_path: Path) -> None:
    with zipfile.ZipFile(archive_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(source_dir.rglob("*")):
            if path.is_dir():
                continue
            arcname = source_dir.name / path.relative_to(source_dir)
            zf.write(path, arcname=str(arcname).replace("\\", "/"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build ESL-PSC toolkit archive plus checksum and manifest."
    )
    parser.add_argument("--dest-dir", required=True, help="Output directory for packaged artifacts")
    parser.add_argument("--rust-main", required=True, help="Path to the esl-psc binary")
    parser.add_argument(
        "--platform",
        required=True,
        choices=("posix", "windows"),
        help="Target wrapper platform family",
    )
    parser.add_argument(
        "--os-tag",
        required=True,
        choices=("linux", "macos", "windows"),
        help="Operating-system tag for artifact naming",
    )
    parser.add_argument("--arch", default="x86_64", help="CPU architecture tag for artifact naming")
    parser.add_argument("--version", help="Release version for artifact naming")
    parser.add_argument("--git-sha", default="", help="Optional source commit SHA for manifest")
    parser.add_argument("--package-name", default="esl-psc-toolkit", help="Artifact base name")
    parser.add_argument(
        "--staging-root",
        default="esl-psc-toolkit",
        help="Top-level directory name inside the archive",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    dest_dir = Path(args.dest_dir).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    version = sanitize_version(args.version) if args.version else sanitize_version(detect_version(repo_root))
    extension = ".zip" if args.platform == "windows" else ".tar.gz"
    archive_name = f"{args.package_name}-v{version}-{args.os_tag}-{args.arch}{extension}"
    archive_path = dest_dir / archive_name

    with tempfile.TemporaryDirectory(prefix="esl-psc-toolkit-") as temp_root:
        stage_dir = Path(temp_root) / args.staging_root
        stage_toolkit(
            dest=stage_dir,
            rust_main=Path(args.rust_main),
            platform=args.platform,
            repo_root=repo_root,
        )
        if args.platform == "windows":
            make_zip(stage_dir, archive_path)
        else:
            make_tar_gz(stage_dir, archive_path)

    checksum = sha256sum(archive_path)
    checksum_path = dest_dir / f"{archive_name}.sha256"
    checksum_path.write_text(f"{checksum}  {archive_name}\n", encoding="utf-8")

    manifest = {
        "package_name": args.package_name,
        "version": version,
        "os": args.os_tag,
        "arch": args.arch,
        "platform_family": args.platform,
        "archive_file": archive_name,
        "archive_sha256": checksum,
        "archive_size_bytes": archive_path.stat().st_size,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_git_sha": args.git_sha,
    }
    manifest_path = dest_dir / f"{archive_name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(str(archive_path))
    print(str(checksum_path))
    print(str(manifest_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
