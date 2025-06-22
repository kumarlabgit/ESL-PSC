"""
Color definitions shared by multiple widgets.
"""

from __future__ import annotations

ZAPPO_STATIC_COLORS: dict[str, str] = {
    # Aliphatic / Hydrophobic
    "I": "#FFC0CB",
    "L": "#FFB5C0",
    "V": "#FFABB6",
    "A": "#FFA0AB",
    "M": "#FF95A0",
    # Aromatic
    "F": "#F1A540",
    "W": "#8B4513",
    "Y": "#FF8F00",
    # Positive
    "K": "#0000CD",
    "R": "#00008B",
    "H": "#4B0082",
    # Negative
    "D": "#FF0000",
    "E": "#E00000",
    # Hydrophilic
    "S": "#90EE90",
    "T": "#86E486",
    "N": "#7CDA7C",
    "Q": "#72D072",
    # Conformationally special
    "P": "#FF00FF",
    "G": "#A000A0",
    # Cysteine
    "C": "#FFFF00",
    # Gap / unknown
    "-": "#C8C8C8",
}
