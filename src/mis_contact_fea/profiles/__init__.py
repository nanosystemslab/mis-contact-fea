"""Profile generators and loaders.

The `SHAPES` dict re-exported here is the canonical registry of built-in
generators (sphere, oblate, prolate, cone, cap, capcone, doublecone).
External profile loaders (CSV, SVG, STEP …) will land in this package
in Phase 4 and be dispatched through `load_profile()`."""
from __future__ import annotations

from .builtins import (
    SHAPES,
    gen_cap_axisym,
    gen_capcone_axisym,
    gen_cone_axisym,
    gen_doublecone_axisym,
    gen_ellipse_axisym,
)


__all__ = [
    "SHAPES",
    "gen_ellipse_axisym",
    "gen_cap_axisym",
    "gen_capcone_axisym",
    "gen_cone_axisym",
    "gen_doublecone_axisym",
]
