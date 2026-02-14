#!/usr/bin/env python3
"""
ESL-PSC GUI Wizard

A graphical interface for running ESL-PSC analyses with an intuitive wizard.
"""
import sys
import os
import signal
from gui.core.logging_utils import setup_logging
import traceback

_NON_WINDOW_QT_PLATFORMS = {"offscreen", "minimal", "minimalegl", "linuxfb", "vnc"}


def _stderr(msg: str) -> None:
    """Write directly to stderr (not routed through GUI debug logging)."""
    sys.__stderr__.write(msg + "\n")
    sys.__stderr__.flush()


def _normalize_qt_platform_env():
    """Prefer a real GUI backend when a display is available."""
    platform = (os.environ.get("QT_QPA_PLATFORM") or "").strip().lower()
    if not platform:
        return

    has_display = bool(os.environ.get("DISPLAY"))
    has_wayland = bool(os.environ.get("WAYLAND_DISPLAY") or os.environ.get("WAYLAND_SOCKET"))
    if has_display or has_wayland:
        if platform in {"offscreen", "minimal", "minimalegl", "linuxfb", "vnc"}:
            os.environ.pop("QT_QPA_PLATFORM", None)
            print(
                f"Cleared QT_QPA_PLATFORM={platform} because a display is available; "
                "using default windowed Qt backend."
            )


def _enforce_windowed_platform(app) -> None:
    """Abort early if Qt selected a non-window backend unexpectedly."""
    if os.environ.get("ESL_PSC_ALLOW_HEADLESS_GUI", "0") == "1":
        return

    platform_name = ""
    screen_count = 0
    try:
        from PySide6.QtGui import QGuiApplication
        platform_name = (QGuiApplication.platformName() or "").strip().lower()
        screen_count = len(QGuiApplication.screens())
    except Exception:
        return

    if platform_name in _NON_WINDOW_QT_PLATFORMS:
        _stderr(
            "Qt selected a non-window platform backend "
            f"('{platform_name}'), so the GUI cannot be shown."
        )
        _stderr(
            "Check your display session and environment variables: "
            "DISPLAY, WAYLAND_DISPLAY, QT_QPA_PLATFORM."
        )
        _stderr(
            "If you intentionally want headless behavior for debugging, set "
            "ESL_PSC_ALLOW_HEADLESS_GUI=1."
        )
        try:
            app.quit()
        except Exception:
            pass
        raise RuntimeError(
            f"Refusing to run GUI on non-window Qt platform '{platform_name}'."
        )

    if screen_count == 0:
        _stderr("Qt started but reports zero screens; GUI cannot be displayed.")
        _stderr(
            "Verify your desktop session (DISPLAY/WAYLAND_DISPLAY) and retry."
        )
        try:
            app.quit()
        except Exception:
            pass
        raise RuntimeError("Qt reports zero screens.")


def main():
    # Delay GUI import until after logging setup to ensure debug prints are redirected
    """Main entry point for the ESL-PSC GUI application."""
    # Initialize application-wide logging and patch print() for GUI modules
    setup_logging()
    _normalize_qt_platform_env()

    from PySide6.QtWidgets import QApplication, QMessageBox
    from PySide6.QtGui import QIcon
    from PySide6.QtCore import Qt, QTimer
    from gui.ui.main_window import MainWindow

    print("Starting ESL-PSC GUI...")
    print(f"Python executable: {sys.executable}")
    print(f"Working directory: {os.getcwd()}")
    
    try:
        # Try to enable high DPI scaling if available without emitting deprecation warnings
        import warnings
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                QApplication.setAttribute(Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
                QApplication.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, True)
                print("Enabled high DPI scaling")
        except (AttributeError, NameError) as e:
            print(f"Could not enable high DPI scaling (not supported in this Qt version): {e}")
        
        # Create application instance
        print("Creating QApplication...")
        app = QApplication(sys.argv)
        _enforce_windowed_platform(app)
        print("QApplication created")

        # ──────────────────────────────────────────────────────────────────
        # Set application icon (affects Dock / task-bar and all dialogs)
        # ──────────────────────────────────────────────────────────────────
        try:
            from pathlib import Path
            icon_path = Path(__file__).resolve().parent.parent / "assets" / "icons" / "app.png"
            if icon_path.is_file():
                app.setWindowIcon(QIcon(str(icon_path)))
                print(f"Set window icon to {icon_path}")
            else:
                print(f"Icon not found at {icon_path}, using default icon")
        except Exception as e:
            print(f"Could not set window icon: {e}")
        
        # Set application metadata
        app.setApplicationName("ESL-PSC Wizard")
        app.setApplicationVersion("2.2.0")
        app.setOrganizationName("ESL-PSC")

        # Make Ctrl+C terminate the Qt event loop when launched from a terminal.
        signal.signal(signal.SIGINT, lambda *_: app.quit())
        sigint_timer = QTimer()
        sigint_timer.timeout.connect(lambda: None)
        sigint_timer.start(200)
        
        # If we are running inside a packaged (Nuitka) build, kick off a
        # background warm-up run of the CLI helper.  This forces the one-file
        # binary to unpack while the user is filling out the wizard so the
        # later real invocation is instantaneous.
        if getattr(sys, "frozen", False) and os.name != 'nt':
            import threading, subprocess
            from pathlib import Path

            def _warm_up_cli():
                try:
                    launcher = Path(os.path.realpath(sys.argv[0]))
                    exe = launcher.with_name("esl_multimatrix" + (".exe" if os.name == "nt" else ""))
                    if exe.is_file():
                        bundle_dir = str(launcher.parent)  # Contents/MacOS or exe dir
                        subprocess.run([
                            str(exe),
                            "--esl_main_dir", bundle_dir,
                            "--help",
                        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                except Exception as e:
                    # Log but never interrupt the GUI
                    print(f"[warm-up] CLI pre-extract failed: {e!r}", file=sys.__stderr__)

            threading.Thread(target=_warm_up_cli, daemon=True).start()

        # Create and show main window
        print("Creating main window...")
        window = MainWindow()
        print("Main window created")
        
        # Show the window
        window.show()
        print("Main window shown")
        
        print("Starting event loop...")
        # Start the event loop
        sys.exit(app.exec())
        
    except Exception as e:
        _stderr(f"FATAL ERROR: {str(e)}")
        traceback.print_exc()
        
        # Try to show error in a message box only when a visible GUI is available.
        try:
            from PySide6.QtGui import QGuiApplication
            platform_name = (QGuiApplication.platformName() or "").strip().lower()
            has_screen = bool(QGuiApplication.primaryScreen())
            if has_screen and platform_name not in _NON_WINDOW_QT_PLATFORMS:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Icon.Critical)
                msg.setText("Fatal Error")
                msg.setInformativeText(str(e))
                msg.setDetailedText(traceback.format_exc())
                msg.setWindowTitle("ESL-PSC GUI Error")
                msg.exec()
        except:
            pass  # If we can't show the error in a GUI, we've already printed it
            
        sys.exit(1)

if __name__ == "__main__":
    main()
