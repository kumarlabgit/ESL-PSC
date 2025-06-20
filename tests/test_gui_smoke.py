# tests/test_gui_smoke.py

from PyQt6.QtWidgets import QApplication

from gui.main import MainWindow

def test_gui_launches_without_error(qt_app):
    """
    Tests if the main GUI window and all its child widgets can be instantiated
    without raising an exception. This is a basic smoke test for the UI.
    
    With the QT_QPA_PLATFORM=offscreen environment variable set, no actual
    window will be shown, but the test will fail if the UI code has a fatal error.
    """
    # 1. The setup is handled by the qt_app fixture.
    
    # 2. The test action is to simply create an instance of our main window.
    # If there are any major errors in the UI construction, this line will fail.
    window = MainWindow()

    # 3. The assertion is that the window object was created successfully.
    assert window is not None
    
    # We can briefly "show" the window to trigger any layout-related code,
    # which is safe to do in the offscreen environment.
    window.show()
    window.close()