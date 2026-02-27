#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import urllib.request
from pathlib import Path


FORMULA_TEMPLATE = """class EslPsc < Formula
  desc "Evolutionary Sparse Learning with Paired Species Contrast"
  homepage "https://github.com/John-Allard/ESL-PSC"
  license "NOASSERTION"
  version "{version}"

  on_macos do
    url "{macos_url}"
    sha256 "{macos_sha}"
  end

  on_linux do
    url "{linux_url}"
    sha256 "{linux_sha}"
  end

  depends_on "python@3.12"
  depends_on "numpy"
  depends_on "pandas"
  depends_on "matplotlib"
  depends_on "seaborn"
  depends_on "biopython"

  def install
    toolkit_root = Dir.children(buildpath).find { |entry| (buildpath/entry).directory? }
    odie "unexpected archive layout: toolkit root directory not found" if toolkit_root.nil?
    libexec.install Dir["#{toolkit_root}/*"]
    (bin/"esl-psc").write_env_script libexec/"bin/esl-psc", {{
      ESL_PSC_PYTHONPATH: "#{libexec}/python",
      ESL_PSC_PYTHON: Formula["python@3.12"].opt_bin/"python3",
    }}
  end

  test do
    assert_match "usage", shell_output("#{bin}/esl-psc --help")
    assert_match "usage", shell_output("#{bin}/esl-psc pairs --help")
    assert_match "usage", shell_output("#{bin}/esl-psc site-counter --help")
    assert_match "usage", shell_output("#{bin}/esl-psc plot --help")
  end
end
"""


def github_api_json(url: str) -> dict:
    headers = {"Accept": "application/vnd.github+json"}
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request) as response:
        return response.read().decode("utf-8")


def parse_sha256_line(raw: str, expected_filename: str) -> str:
    line = raw.strip().splitlines()[0].strip()
    parts = line.split()
    if len(parts) < 2:
        raise ValueError(f"invalid sha256 file content: {line!r}")
    checksum = parts[0]
    filename = parts[-1].lstrip("*")
    if filename != expected_filename:
        raise ValueError(
            f"sha256 file mismatch: expected {expected_filename}, found {filename}"
        )
    return checksum


def main() -> int:
    parser = argparse.ArgumentParser(description="Render Homebrew formula from GitHub release assets.")
    parser.add_argument("--repo", required=True, help="GitHub repo in owner/name form")
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v2.4.1")
    parser.add_argument("--output", required=True, help="Output .rb formula path")
    args = parser.parse_args()

    release = github_api_json(f"https://api.github.com/repos/{args.repo}/releases/tags/{args.tag}")
    assets = {asset["name"]: asset["browser_download_url"] for asset in release.get("assets", [])}

    version_match = re.fullmatch(r"v(.+)", args.tag.strip())
    if not version_match:
        raise ValueError(f"invalid tag format: {args.tag}")
    version = version_match.group(1)

    linux_name = f"esl-psc-toolkit-v{version}-linux-x86_64.tar.gz"
    linux_sha_name = f"{linux_name}.sha256"
    macos_name = f"esl-psc-toolkit-v{version}-macos-x86_64.tar.gz"
    macos_sha_name = f"{macos_name}.sha256"

    missing = [name for name in (linux_name, linux_sha_name, macos_name, macos_sha_name) if name not in assets]
    if missing:
        raise KeyError(f"missing release assets for formula generation: {', '.join(missing)}")

    linux_sha = parse_sha256_line(fetch_text(assets[linux_sha_name]), linux_name)
    macos_sha = parse_sha256_line(fetch_text(assets[macos_sha_name]), macos_name)

    rendered = FORMULA_TEMPLATE.format(
        version=version,
        linux_url=assets[linux_name],
        linux_sha=linux_sha,
        macos_url=assets[macos_name],
        macos_sha=macos_sha,
    )

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
