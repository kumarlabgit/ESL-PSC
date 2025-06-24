"""Utility module for application-wide logging configuration.

This centralises logging setup so that all GUI-related debug messages can
be silenced in production builds yet re-enabled easily when troubleshooting.

How it works
============
1. ``setup_logging`` configures the *root* logger with a sane format.
   The log level is DEBUG when the environment variable
   ``ESL_PSC_GUI_DEBUG=1`` is present, otherwise INFO (you can customise
   this by passing ``debug=True/False`` explicitly).

2. Legacy ``print`` statements that live inside *our* GUI package (i.e.
   any module whose ``__name__`` starts with ``"gui"``) are preserved but
   transparently rerouted to the logging system as ``logger.debug`` calls.
   This means:
   • No code changes are required – we keep the existing prints for
     convenience during development.
   • In production, because the log level defaults to INFO, these debug
     prints are effectively hidden.
   • Setting ``ESL_PSC_GUI_DEBUG=1`` (or passing ``debug=True``) brings
     them back instantly.

3. ``print`` calls originating from *outside* the GUI package are left
   untouched.  This is important so that the CLI utilities executed in
   the worker thread (e.g. ``esl_multimatrix``) continue to be captured
   by the existing ``StreamEmitter`` redirection logic and shown in the
   GUI’s terminal pane.
"""
from __future__ import annotations

import builtins
import inspect
import logging
import os
from typing import Any

# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def setup_logging(debug: bool | None = None, *, app_name: str = "ESL-PSC GUI") -> logging.Logger:
    """Initialise global logging.

    Parameters
    ----------
    debug: bool | None
        When *True*, force DEBUG level output.  When *False*, use INFO.
        When *None* (default), the level is decided by the environment
        variable ``ESL_PSC_GUI_DEBUG`` ("1" == debug on).
    app_name: str
        Name of the top-level logger that will receive redirected prints.

    Returns
    -------
    logging.Logger
        The application logger instance (useful for immediate logging).
    """
    if debug is None:
        debug = os.getenv("ESL_PSC_GUI_DEBUG", "0") == "1"

    level = logging.DEBUG if debug else logging.INFO

    # Root logger configuration (only if the user hasn’t configured it yet)
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%H:%M:%S",
        )
    else:
        # Even if handlers exist, still honour the requested level
        logging.getLogger().setLevel(level)

    # Silence extremely verbose third-party modules when debugging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.INFO)

    app_logger = logging.getLogger(app_name)

    _patch_print(app_logger)
    return app_logger

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _patch_print(app_logger: logging.Logger) -> None:
    """Monkey-patch :pyfunc:`print` so that GUI debug messages go via logging.

    Only modules whose ``__name__`` begins with ``"gui"`` are redirected.
    All other calls fall back to the original builtin and therefore behave
    exactly as before (allowing the existing worker stream interception to
    keep functioning).
    """
    orig_print = builtins.print  # Preserve original for non-GUI callers

    # Avoid double-patching
    if getattr(builtins, "_esl_gui_print_patched", False):
        return

    def _print_proxy(*args: Any, **kwargs: Any) -> None:  # noqa: D401 – simple proxy
        # Cheaply determine caller module name
        frame = inspect.currentframe()
        if frame is None:
            orig_print(*args, **kwargs)
            return
        caller_frame = frame.f_back
        module_name = caller_frame.f_globals.get("__name__", "") if caller_frame else ""

        is_gui_module = module_name.startswith("gui")
        # Special-case when the entry point is executed as "-m gui.main" where
        # the module name becomes "__main__". We look at the caller's file path
        # to decide if it lives inside our `gui` package so that those prints
        # are still considered GUI debug messages (and therefore can be hidden
        # in production).
        if not is_gui_module and module_name == "__main__":
            caller_file = caller_frame.f_globals.get("__file__", "")
            try:
                # Normalise path for robust substring check
                import os as _os
                caller_file_norm = _os.path.abspath(caller_file)
                is_gui_module = _os.sep + "gui" + _os.sep in caller_file_norm
            except Exception:
                # If any error occurs during path inspection, fall back to False
                is_gui_module = False

        if is_gui_module:
            # Build message similarly to the builtin print (respect sep/end)
            sep: str = kwargs.get("sep", " ")
            end: str = kwargs.get("end", "\n")
            message = sep.join(str(a) for a in args) + end.rstrip("\n")
            app_logger.debug(message)
        else:
            # Not our code – delegate untouched
            orig_print(*args, **kwargs)

    builtins.print = _print_proxy  # type: ignore[assignment]
    builtins._esl_gui_print_patched = True
