#!/usr/bin/env python3
import argparse
import re
import subprocess
from pathlib import Path


def collect_python_flags(repo_root: Path):
    import importlib.util
    import sys

    sys.path.insert(0, str(repo_root))
    from esl_psc_cli import esl_integrator as ei
    from esl_psc_cli import deletion_canceler as dc

    parser = argparse.ArgumentParser()
    ei.get_esl_args(parser)
    ei.get_single_run_esl_args(parser, input_alignments_req=False)
    dc.get_deletion_canceler_args(parser)

    # multimatrix-local args are defined in esl_multimatrix.main()
    # and not exposed via helper; parse source for add_argument calls.
    mm_src = (repo_root / "esl_psc_cli" / "esl_multimatrix.py").read_text()
    mm_flags = set(re.findall(r"add_argument\(\s*'(--[a-zA-Z0-9_\-]+)'", mm_src))

    flags = set(mm_flags)
    for action in parser._actions:
        for opt in action.option_strings:
            if opt.startswith("--"):
                flags.add(opt)

    flags.discard("--help")
    return flags


def collect_rust_flags(repo_root: Path, bin_path: str | None):
    if bin_path:
        cmd = [bin_path, "--help"]
        cwd = repo_root
    else:
        cmd = ["cargo", "run", "--quiet", "--manifest-path", str(repo_root / "esl_unified_rs" / "Cargo.toml"), "--", "--help"]
        cwd = repo_root

    out = subprocess.check_output(cmd, cwd=cwd, text=True)
    flags = set(re.findall(r"--[a-zA-Z0-9_\-]+", out))
    flags.discard("--help")
    return flags


def main():
    p = argparse.ArgumentParser(description="Check CLI flag parity between Python ESL-PSC and Rust unified CLI")
    p.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    p.add_argument("--rust-bin", default=None, help="Path to esl_unified_rs binary (optional)")
    args = p.parse_args()

    root = Path(args.repo_root).resolve()
    py_flags = collect_python_flags(root)
    rs_flags = collect_rust_flags(root, args.rust_bin)

    missing = sorted(py_flags - rs_flags)
    extra = sorted(rs_flags - py_flags)

    print(f"python_flags={len(py_flags)} rust_flags={len(rs_flags)}")
    print(f"missing_in_rust={len(missing)}")
    for f in missing:
        print(f"  MISSING {f}")
    print(f"extra_in_rust={len(extra)}")
    for f in extra:
        print(f"  EXTRA {f}")

    if missing:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
