import math

from gui.ui.widgets.graphics_view import ZoomableGraphicsView


def test_trackpad_zoom_step_is_damped(qt_app):
    view = ZoomableGraphicsView()

    # A typical high-resolution trackpad delta should produce a small zoom
    # change, not a full mouse-wheel notch.
    steps = 120 / view._PIXELS_PER_TRACKPAD_STEP
    capped_steps = min(view._MAX_TRACKPAD_STEPS_PER_EVENT, steps)
    view._apply_zoom_steps(capped_steps)

    assert view._zoom_level <= view._MAX_TRACKPAD_STEPS_PER_EVENT
    assert math.isclose(view.transform().m11(), math.pow(1.15, capped_steps), rel_tol=1e-6)
    view.close()


def test_mouse_wheel_zoom_step_remains_one_notch(qt_app):
    view = ZoomableGraphicsView()

    view._apply_zoom_steps(1)

    assert view._zoom_level == 1
    assert math.isclose(view.transform().m11(), 1.15, rel_tol=1e-6)
    view.close()

