from __future__ import annotations

from typing import List, Tuple
from PySide6.QtGui import QColor


def viridis_stops_rgb() -> List[Tuple[float, Tuple[int, int, int]]]:
    """Return Viridis control points as (t, (r,g,b)) tuples.

    Stops are approximate Matplotlib Viridis stops and are suitable for
    building QLinearGradient or for interpolation.
    """
    return [
        (0.0,  (68, 1, 84)),     # #440154
        (0.13, (72, 40, 120)),   # #482878
        (0.25, (62, 73, 137)),   # #3e4989
        (0.38, (49, 104, 142)),  # #31688e
        (0.50, (38, 130, 142)),  # #26828e
        (0.63, (31, 158, 137)),  # #1f9e89
        (0.75, (53, 183, 121)),  # #35b779
        (0.88, (143, 215, 68)),  # #8fd744
        (1.0,  (253, 231, 37)),  # #fde725
    ]


def viridis_qcolor(p: float) -> QColor:
    """Interpolate a QColor from the Viridis colormap for p in [0, 1]."""
    if p < 0.0:
        p = 0.0
    elif p > 1.0:
        p = 1.0
    stops = viridis_stops_rgb()
    if p <= stops[0][0]:
        r, g, b = stops[0][1]
        return QColor(r, g, b)
    if p >= stops[-1][0]:
        r, g, b = stops[-1][1]
        return QColor(r, g, b)
    for i in range(len(stops) - 1):
        t1, c1 = stops[i]
        t2, c2 = stops[i + 1]
        if t1 <= p <= t2:
            if t2 == t1:
                r, g, b = c2
                return QColor(r, g, b)
            a = (p - t1) / (t2 - t1)
            r = int(round(c1[0] + a * (c2[0] - c1[0])))
            g = int(round(c1[1] + a * (c2[1] - c1[1])))
            b = int(round(c1[2] + a * (c2[2] - c1[2])))
            return QColor(r, g, b)
    # Fallback (shouldn't reach)
    r, g, b = stops[-1][1]
    return QColor(r, g, b)
