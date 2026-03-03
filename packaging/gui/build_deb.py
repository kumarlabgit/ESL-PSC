#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]

CONTROL_TEMPLATE = """Package: esl-psc-gui
Version: {version}
Section: science
Priority: optional
Architecture: {arch}
Maintainer: ESL-PSC Maintainers <maintainers@esl-psc.local>
Description: ESL-PSC graphical interface
 Evolutionary Sparse Learning with Paired Species Contrast (ESL-PSC)
 graphical desktop application with bundled runtime helpers.
"""

DESKTOP_ENTRY = """[Desktop Entry]
Type=Application
Name=ESL-PSC
Exec=esl-psc-gui
Icon=esl-psc
Categories=Science;
Terminal=false
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


def find_gui_executable_name(dist_dir: Path) -> str:
    preferred_names = ("main.bin", "main", "ESL-PSC", "ESL-PSC.bin")
    for name in preferred_names:
        candidate = dist_dir / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return name

    for candidate in dist_dir.iterdir():
        if not candidate.is_file():
            continue
        if os.access(candidate, os.X_OK):
            return candidate.name

    raise FileNotFoundError(f"unable to find executable in {dist_dir}")


def write_sha256(path: Path) -> Path:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    sha_path = Path(str(path) + ".sha256")
    sha_path.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return sha_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Debian package for ESL-PSC GUI.")
    parser.add_argument("--dist-dir", required=True, help="Path to Nuitka .dist payload directory")
    parser.add_argument("--dest-dir", required=True, help="Output directory for .deb artifact")
    parser.add_argument("--version", help="Package version (default: from esl_psc_rs/Cargo.toml)")
    parser.add_argument("--arch", default="amd64", help="Debian architecture (amd64|arm64)")
    args = parser.parse_args()

    dist_dir = Path(args.dist_dir).resolve()
    if not dist_dir.is_dir():
        raise FileNotFoundError(f"dist directory not found: {dist_dir}")

    repo_root = REPO_ROOT
    dest_dir = Path(args.dest_dir).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    version = args.version or detect_version(repo_root)
    arch = normalize_arch(args.arch)
    exe_name = find_gui_executable_name(dist_dir)

    with tempfile.TemporaryDirectory(prefix="esl-psc-gui-deb-") as tmp_root:
        tmp_root_path = Path(tmp_root)
        pkg_root = tmp_root_path / f"esl-psc-gui_{version}_{arch}"
        debian_dir = pkg_root / "DEBIAN"
        opt_dir = pkg_root / "opt" / "esl-psc-gui"
        usr_bin = pkg_root / "usr" / "bin"
        app_dir = pkg_root / "usr" / "share" / "applications"
        icon_dir = pkg_root / "usr" / "share" / "icons" / "hicolor" / "256x256" / "apps"
        doc_dir = pkg_root / "usr" / "share" / "doc" / "esl-psc-gui"

        debian_dir.mkdir(parents=True, exist_ok=True)
        opt_dir.mkdir(parents=True, exist_ok=True)
        usr_bin.mkdir(parents=True, exist_ok=True)
        app_dir.mkdir(parents=True, exist_ok=True)
        icon_dir.mkdir(parents=True, exist_ok=True)
        doc_dir.mkdir(parents=True, exist_ok=True)

        shutil.copytree(dist_dir, opt_dir, dirs_exist_ok=True)

        launcher_path = usr_bin / "esl-psc-gui"
        launcher_path.write_text(
            "\n".join(
                [
                    "#!/bin/sh",
                    f'exec "/opt/esl-psc-gui/{exe_name}" "$@"',
                    "",
                ]
            ),
            encoding="utf-8",
        )
        launcher_path.chmod(0o755)

        (app_dir / "esl-psc.desktop").write_text(DESKTOP_ENTRY, encoding="utf-8")
        shutil.copy2(repo_root / "assets" / "icons" / "app.png", icon_dir / "esl-psc.png")
        shutil.copy2(repo_root / "gui" / "README.md", doc_dir / "README.md")

        control_text = CONTROL_TEMPLATE.format(version=version, arch=arch)
        (debian_dir / "control").write_text(control_text, encoding="utf-8")

        out_deb = dest_dir / f"esl-psc-gui_{version}_{arch}.deb"
        run_checked(["dpkg-deb", "--build", str(pkg_root), str(out_deb)])
        sha_path = write_sha256(out_deb)
        print(out_deb)
        print(sha_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
