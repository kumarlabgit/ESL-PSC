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

# pytest-qt provides a `qapp` fixture; in minimal environments it may be absent.
# Our GUI tests only need a QApplication instance + processEvents().
try:
    import pytestqt  # type: ignore  # noqa: F401
    _HAVE_PYTEST_QT = True
except Exception:
    _HAVE_PYTEST_QT = False

if not _HAVE_PYTEST_QT:
    @pytest.fixture(scope="session")
    def qapp():  # type: ignore
        from PySide6.QtWidgets import QApplication
        app = QApplication.instance()
        if app is None:
            app = QApplication([])
        return app

# pytest-qt 4 uses `qapp`; older code may expect `qt_app`.
@pytest.fixture
def qt_app(qapp):  # type: ignore
    """Alias for the qapp fixture to keep legacy tests working."""
    return qapp
