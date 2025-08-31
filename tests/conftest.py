"""pytest configuration to ensure project root is on sys.path.

Allows `import gui` and other top-level packages during test discovery.
"""
import os
import sys
from pathlib import Path

# Force Qt to use an offscreen platform so GUI smoke tests run in headless CI
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

import pytest

# pytest-qt 4 uses `qapp`; older code may expect `qt_app`.
@pytest.fixture
def qt_app(qapp):  # type: ignore
    """Alias for the qapp fixture to keep legacy tests working."""
    return qapp
