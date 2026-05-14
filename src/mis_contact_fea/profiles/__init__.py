"""Profile generators and loaders.

The `SHAPES` dict is the canonical registry of built-in generators
(sphere, oblate, prolate, cone, cap, capcone, doublecone). External
loaders for CSV / SVG / STEP files live alongside in `from_csv`,
`from_svg`, `from_step` and are dispatched through `load_profile()`.

Internal canonical unit is **micrometre**. Loaders convert from their
declared input unit at load time, so the rest of the pipeline always
sees µm.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from .builtins import (
    SHAPES,
    gen_cap_axisym,
    gen_capcone_axisym,
    gen_cone_axisym,
    gen_doublecone_axisym,
    gen_ellipse_axisym,
)


if TYPE_CHECKING:  # avoid a circular import at runtime
    from mis_contact_fea.config import ProfileConfig


# Length-unit conversion: factor applied to convert input -> µm.
UNIT_TO_UM: dict[str, float] = {
    "um": 1.0,
    "mm": 1_000.0,
    "cm": 10_000.0,
    "m":  1_000_000.0,
}


def um_per_unit(unit: str) -> float:
    """Return how many micrometres are in one of the given unit."""
    try:
        return UNIT_TO_UM[unit]
    except KeyError as e:
        raise ValueError(
            f"unknown unit {unit!r}; must be one of {list(UNIT_TO_UM)}"
        ) from e


def load_profile(
    spec: "ProfileConfig",
    Hp: float = 220.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the right-half (x, z) polyline for the configured profile.

    Dispatches on `spec.kind`:

    - `builtin` — look up `spec.name` in `SHAPES`
    - `csv`     — `from_csv.load_csv_profile(spec.path, units=spec.units, ...)`
    - `svg`     — `from_svg.load_svg_profile(spec.path, units=spec.units, ...)`
    - `step`    — `from_step.load_step_profile(spec.path, units=spec.units, ...)`

    All units are converted to µm before return."""
    kind = spec.kind
    if kind == "builtin":
        if not spec.name or spec.name not in SHAPES:
            raise ValueError(
                f"profile.kind=builtin requires a `name` in {list(SHAPES)}; got {spec.name!r}"
            )
        return SHAPES[spec.name](Hp=Hp)

    if spec.path is None:
        raise ValueError(f"profile.kind={kind!r} requires a `path` field")

    if kind == "csv":
        from .from_csv import load_csv_profile
        return load_csv_profile(
            spec.path, Hp=Hp, units=spec.units,
            prepend_pillar=spec.prepend_pillar,
        )
    if kind == "svg":
        from .from_svg import load_svg_profile
        return load_svg_profile(
            spec.path, Hp=Hp, units=spec.units,
            prepend_pillar=spec.prepend_pillar,
        )
    if kind == "step":
        from .from_step import load_step_profile
        return load_step_profile(spec.path, Hp=Hp, units=spec.units)

    raise ValueError(f"unknown profile kind: {kind!r}")


__all__ = [
    "SHAPES",
    "UNIT_TO_UM",
    "um_per_unit",
    "load_profile",
    "gen_ellipse_axisym",
    "gen_cap_axisym",
    "gen_capcone_axisym",
    "gen_cone_axisym",
    "gen_doublecone_axisym",
]
