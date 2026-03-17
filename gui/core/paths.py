from __future__ import annotations

import os
from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def artifacts_dir() -> Path:
    return project_root() / ".artifacts"


def artifact_path(*parts: str) -> str:
    return str(artifacts_dir().joinpath(*parts))


def default_output_dir() -> str:
    """Return a cross-platform default for user-visible analysis outputs."""
    documents_dir = ""
    home_dir = os.path.expanduser("~")
    try:
        from PySide6.QtCore import QStandardPaths

        documents_dir = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DocumentsLocation)
        if not documents_dir:
            documents_dir = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.HomeLocation)
    except Exception:
        documents_dir = ""

    base_dir = documents_dir or home_dir
    return str(Path(base_dir) / "ESL-PSC")
