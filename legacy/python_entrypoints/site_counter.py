#!/usr/bin/env python3
import sys

from esl_psc_cli.fast_scan_cli import main


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
