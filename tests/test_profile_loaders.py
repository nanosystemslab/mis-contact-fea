"""Tests for the profile loaders and the load_profile() dispatch.

Round-trip a built-in shape through the CSV loader at every supported
unit, confirm the dispatch routes correctly, and check that the
chatter-prone "must start/end on axis" validation fires for malformed
input.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from mis_contact_fea.config import ProfileConfig
from mis_contact_fea.profiles import SHAPES, UNIT_TO_UM, load_profile, um_per_unit
from mis_contact_fea.profiles.from_csv import load_csv_profile


def _write_sphere_as_csv(tmp_path: Path, unit: str) -> Path:
    """Write the sphere's (x, z) polyline as a CSV in the given unit."""
    x_um, z_um = SHAPES["sphere"](Hp=220.0)
    scale = UNIT_TO_UM[unit]  # µm per unit; we DIVIDE by this when writing
    csv_path = tmp_path / f"sphere_in_{unit}.csv"
    with csv_path.open("w") as f:
        w = csv.writer(f)
        w.writerow([f"x_{unit}", f"z_{unit}"])
        for xi, zi in zip(x_um / scale, z_um / scale):
            w.writerow([f"{xi:.10g}", f"{zi:.10g}"])
    return csv_path


@pytest.mark.parametrize("unit", ["um", "mm", "cm", "m"])
def test_csv_roundtrip_all_units(tmp_path: Path, unit: str) -> None:
    """Save sphere as CSV in `unit`, load via dispatch, compare to original."""
    x_orig, z_orig = SHAPES["sphere"](Hp=220.0)
    csv_path = _write_sphere_as_csv(tmp_path, unit)

    spec = ProfileConfig(kind="csv", path=csv_path, units=unit)
    x_loaded, z_loaded = load_profile(spec, Hp=220.0)

    assert len(x_loaded) == len(x_orig), "polyline length mismatch"
    np.testing.assert_allclose(x_loaded, x_orig, atol=1e-3,
                               err_msg=f"x mismatch for unit={unit}")
    np.testing.assert_allclose(z_loaded, z_orig, atol=1e-3,
                               err_msg=f"z mismatch for unit={unit}")


def test_builtin_dispatch() -> None:
    """load_profile() with kind=builtin should equal SHAPES[name]() directly."""
    spec = ProfileConfig(kind="builtin", name="cone")
    x_d, z_d = load_profile(spec, Hp=220.0)
    x_b, z_b = SHAPES["cone"](Hp=220.0)
    np.testing.assert_array_equal(x_d, x_b)
    np.testing.assert_array_equal(z_d, z_b)


def test_unknown_unit_rejected() -> None:
    with pytest.raises(ValueError, match="unknown unit"):
        um_per_unit("furlongs")


def test_csv_axis_endpoint_validation(tmp_path: Path) -> None:
    """Polyline that doesn't start on the axis should be rejected."""
    p = tmp_path / "bad.csv"
    p.write_text("x_um,z_um\n10,0\n20,10\n0,20\n")
    with pytest.raises(ValueError, match="first point must be on the axis"):
        load_csv_profile(p)


def test_csv_negative_x_rejected(tmp_path: Path) -> None:
    """Right-half-only constraint: x must be >= 0."""
    p = tmp_path / "neg.csv"
    p.write_text("x_um,z_um\n0,-220\n-5,0\n0,50\n")
    with pytest.raises(ValueError, match="x must be >= 0"):
        load_csv_profile(p)
