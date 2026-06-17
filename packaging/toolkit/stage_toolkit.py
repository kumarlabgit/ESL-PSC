#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def copy_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def remove_bytecode(root: Path) -> None:
    for pycache_dir in root.rglob("__pycache__"):
        if pycache_dir.is_dir():
            shutil.rmtree(pycache_dir)
    for pyc in root.rglob("*.pyc"):
        pyc.unlink(missing_ok=True)


def stage_toolkit(
    *,
    dest: Path,
    rust_main: Path,
    platform: str,
    repo_root: Path | None = None,
) -> None:
    if platform not in ("posix", "windows"):
        raise ValueError(f"unsupported platform family: {platform}")
    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[2]

    dest = dest.resolve()
    rust_main = rust_main.resolve()
    bin_dir = dest / "bin"
    py_root = dest / "python"
    gui_core_dir = py_root / "gui" / "core"

    bin_dir.mkdir(parents=True, exist_ok=True)
    gui_core_dir.mkdir(parents=True, exist_ok=True)

    if not rust_main.is_file():
        raise FileNotFoundError(f"missing Rust main binary: {rust_main}")

    shutil.copy2(rust_main, bin_dir / rust_main.name)

    copy_tree(repo_root / "esl_psc_cli", py_root / "esl_psc_cli")

    shutil.copy2(repo_root / "gui" / "__init__.py", py_root / "gui" / "__init__.py")
    shutil.copy2(repo_root / "gui" / "core" / "fast_scan.py", gui_core_dir / "fast_scan.py")
    shutil.copy2(repo_root / "gui" / "core" / "fasta_io.py", gui_core_dir / "fasta_io.py")
    shutil.copy2(
        repo_root / "gui" / "core" / "ancestral_reconstruction.py",
        gui_core_dir / "ancestral_reconstruction.py",
    )

    shutil.copy2(repo_root / "packaging" / "toolkit" / "README.md", dest / "README.md")
    shutil.copy2(repo_root / "requirements-toolkit.txt", dest / "requirements-toolkit.txt")

    remove_bytecode(py_root)


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage ESL-PSC toolkit payload.")
    parser.add_argument("--dest", required=True, help="Toolkit output directory")
    parser.add_argument("--rust-main", required=True, help="Path to esl-psc binary")
    parser.add_argument(
        "--platform",
        required=True,
        choices=("posix", "windows"),
        help="Target platform family for wrapper scripts",
    )
    args = parser.parse_args()
    stage_toolkit(
        dest=Path(args.dest),
        rust_main=Path(args.rust_main),
        platform=args.platform,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
