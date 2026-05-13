"""Smoke tests for mis_contact_fea.profiles — verify every SHAPES entry generates a
non-empty polyline and that the effective outer diameter is uniform to within
10 nm across all shapes (the cone/doublecone fillet calibration tolerance)."""
from __future__ import annotations

import numpy as np

from mis_contact_fea.profiles import SHAPES


TARGET_DIAMETER = 100.0  # µm
TOLERANCE = 0.01  # µm — 10 nm; fillet calibration is to ~4 nm


def test_all_shapes_generate_polylines() -> None:
    for name, fn in SHAPES.items():
        x, z = fn()
        assert len(x) > 5, f"{name}: polyline too short"
        assert len(x) == len(z), f"{name}: x/z length mismatch"
        assert np.isfinite(x).all() and np.isfinite(z).all(), f"{name}: non-finite point"


def test_uniform_effective_diameter() -> None:
    diameters = {}
    for name, fn in SHAPES.items():
        x, _ = fn()
        diameters[name] = 2 * float(x.max())

    for name, d in diameters.items():
        assert abs(d - TARGET_DIAMETER) < TOLERANCE, (
            f"{name}: effective diameter {d:.4f} µm differs from target "
            f"{TARGET_DIAMETER} by more than {TOLERANCE} µm"
        )


def test_polyline_starts_on_axis() -> None:
    for name, fn in SHAPES.items():
        x, _ = fn()
        assert abs(x[0]) < 1e-9, f"{name}: polyline does not start on x=0 axis"
