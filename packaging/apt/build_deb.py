#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "packaging" / "toolkit"))
from stage_toolkit import stage_toolkit  # type: ignore


CONTROL_TEMPLATE = """Package: esl-psc
Version: {version}
Section: science
Priority: optional
Architecture: {arch}
Depends: python3 (>= 3.10), python3-numpy, python3-pandas, python3-matplotlib, python3-seaborn, python3-biopython
Maintainer: ESL-PSC Maintainers <maintainers@esl-psc.local>
Description: ESL-PSC unified Rust CLI toolkit
 Evolutionary Sparse Learning with Paired Species Contrast (ESL-PSC)
 unified command-line toolkit. Provides:
  - esl-psc run
  - esl-psc pairs
  - esl-psc site-counter
  - esl-psc plot
"""


def run_checked(command: list[str], cwd: Path | None = None) -> None:
    result = subprocess.run(command, cwd=cwd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}")


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


def normalize_arch(arch: str) -> str:
    mapping = {
        "x86_64": "amd64",
        "amd64": "amd64",
        "aarch64": "arm64",
        "arm64": "arm64",
    }
    if arch not in mapping:
        raise ValueError(f"unsupported deb architecture: {arch}")
    return mapping[arch]


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Debian package for esl-psc toolkit.")
    parser.add_argument("--dest-dir", required=True, help="Output directory for .deb artifact")
    parser.add_argument("--rust-main", required=True, help="Path to built esl-psc binary")
    parser.add_argument("--version", help="Package version (default: from esl_psc_rs/Cargo.toml)")
    parser.add_argument("--arch", default="amd64", help="Debian architecture (amd64|arm64)")
    args = parser.parse_args()

    repo_root = REPO_ROOT
    dest_dir = Path(args.dest_dir).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    version = args.version or detect_version(repo_root)
    arch = normalize_arch(args.arch)

    with tempfile.TemporaryDirectory(prefix="esl-psc-deb-") as tmp_root:
        tmp_root_path = Path(tmp_root)
        toolkit_stage = tmp_root_path / "toolkit_stage"
        stage_toolkit(
            dest=toolkit_stage,
            rust_main=Path(args.rust_main),
            platform="posix",
            repo_root=repo_root,
        )

        pkg_root = tmp_root_path / f"esl-psc_{version}_{arch}"
        debian_dir = pkg_root / "DEBIAN"
        usr_bin = pkg_root / "usr" / "bin"
        usr_lib = pkg_root / "usr" / "lib" / "esl-psc"
        usr_share_doc = pkg_root / "usr" / "share" / "doc" / "esl-psc"

        debian_dir.mkdir(parents=True, exist_ok=True)
        usr_bin.mkdir(parents=True, exist_ok=True)
        usr_lib.mkdir(parents=True, exist_ok=True)
        usr_share_doc.mkdir(parents=True, exist_ok=True)

        shutil.copy2(toolkit_stage / "bin" / "esl-psc", usr_bin / "esl-psc")
        (usr_bin / "esl-psc").chmod(0o755)

        shutil.copytree(toolkit_stage / "python", usr_lib / "python", dirs_exist_ok=True)
        shutil.copy2(repo_root / "README.md", usr_share_doc / "README.md")

        control_text = CONTROL_TEMPLATE.format(version=version, arch=arch)
        (debian_dir / "control").write_text(control_text, encoding="utf-8")

        out_deb = dest_dir / f"esl-psc_{version}_{arch}.deb"
        run_checked(["dpkg-deb", "--build", str(pkg_root), str(out_deb)])
        print(out_deb)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
