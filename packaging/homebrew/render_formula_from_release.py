#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import urllib.request
from hashlib import sha256
from pathlib import Path


FORMULA_TEMPLATE = """class EslPsc < Formula
  desc "Evolutionary Sparse Learning with Paired Species Contrast"
  homepage "https://github.com/John-Allard/ESL-PSC"
  license "NOASSERTION"
  url "{source_url}"
  sha256 "{source_sha}"
  version "{version}"

  depends_on "rust" => :build
  depends_on "python@3.12"
  depends_on "numpy"
  depends_on "pandas"
  depends_on "matplotlib"
  depends_on "seaborn"
  depends_on "biopython"

  def install
    system "cargo", "install", "--path", "esl_psc_rs", "--root", libexec

    (libexec/"python/esl_psc_cli").install Dir["esl_psc_cli/*"]
    (libexec/"python/gui").install "gui/__init__.py"
    (libexec/"python/gui/core").mkpath
    (libexec/"python/gui/core").install "gui/core/fast_scan.py"
    (libexec/"python/gui/core").install "gui/core/fasta_io.py"
    (libexec/"python/gui/core").install "gui/core/ancestral_reconstruction.py"

    (bin/"esl-psc").write_env_script libexec/"bin/esl-psc", {{
      ESL_PSC_PYTHONPATH: "#{{libexec}}/python",
      ESL_PSC_PYTHON: Formula["python@3.12"].opt_bin/"python3",
    }}
  end

  test do
    assert_match "usage", shell_output("#{{bin}}/esl-psc --help")
    assert_match "usage", shell_output("#{{bin}}/esl-psc pairs --help")
    assert_match "usage", shell_output("#{{bin}}/esl-psc site-counter --help")
    assert_match "usage", shell_output("#{{bin}}/esl-psc plot --help")
  end
end
"""


def sha256_of_url(url: str) -> str:
    digest = sha256()
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request) as response:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Render Homebrew formula from GitHub source tag tarball.")
    parser.add_argument("--repo", required=True, help="GitHub repo in owner/name form")
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v2.4.1")
    parser.add_argument("--output", required=True, help="Output .rb formula path")
    args = parser.parse_args()

    version_match = re.fullmatch(r"v(.+)", args.tag.strip())
    if not version_match:
        raise ValueError(f"invalid tag format: {args.tag}")
    version = version_match.group(1)

    source_url = f"https://github.com/{args.repo}/archive/refs/tags/{args.tag}.tar.gz"
    source_sha = sha256_of_url(source_url)

    rendered = FORMULA_TEMPLATE.format(
        version=version,
        source_url=source_url,
        source_sha=source_sha,
    )

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
