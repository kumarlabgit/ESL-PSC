#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import re
import urllib.request
from pathlib import Path


def sha256_url(url: str) -> str:
    digest = hashlib.sha256()
    with urllib.request.urlopen(url) as response:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Render conda-forge meta.yaml for a release tag.")
    parser.add_argument("--repo", required=True, help="GitHub repo in owner/name form")
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v2.4.1")
    parser.add_argument(
        "--template",
        default="packaging/conda-forge/recipe/meta.yaml",
        help="Template meta.yaml path",
    )
    parser.add_argument("--output", required=True, help="Output meta.yaml path")
    args = parser.parse_args()

    tag = args.tag.strip()
    match = re.fullmatch(r"v(.+)", tag)
    if not match:
        raise ValueError(f"invalid tag format: {tag}")
    version = match.group(1)

    source_url = f"https://github.com/{args.repo}/archive/refs/tags/{tag}.tar.gz"
    source_sha = sha256_url(source_url)

    template_path = Path(args.template).resolve()
    rendered = template_path.read_text(encoding="utf-8")
    rendered = re.sub(r'(\{% set version = ")[^"]+("\s*%\})', rf"\g<1>{version}\2", rendered)
    rendered = re.sub(
        r'(\{% set sha256 = ")[^"]+("\s*%\})',
        rf"\g<1>{source_sha}\2",
        rendered,
    )

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
