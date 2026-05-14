"""Tests for asymmetric top/bottom profile support.

Covers two regimes:

1. Backward compat — when `top_profile` is omitted (or `bottom_z_apex`
   isn't passed to `top_body_vertices`), behavior is byte-identical to
   the previous single-profile code path.
2. Asymmetric — when the top body uses a different profile, its apex
   still lands exactly `initial_gap` above the bottom apex.
"""
from __future__ import annotations

import numpy as np
import pytest

from mis_contact_fea.config import SimulationConfig
from mis_contact_fea.mesh.layout import bottom_body_vertices, top_body_vertices
from mis_contact_fea.profiles import SHAPES


P = 198.0
PLATE = 30.0
GAP = -55.0


def test_top_body_back_compat_symmetric() -> None:
    """No `bottom_z_apex` kwarg = identical output to legacy call."""
    half_x, half_z = SHAPES["sphere"](Hp=220.0)
    legacy_args = (half_x, half_z, P, PLATE, GAP)
    vx_old, vz_old, tags_old = top_body_vertices(*legacy_args)
    vx_new, vz_new, tags_new = top_body_vertices(*legacy_args, bottom_z_apex=None)
    np.testing.assert_array_equal(vx_old, vx_new)
    np.testing.assert_array_equal(vz_old, vz_new)
    assert tags_old == tags_new


def test_top_apex_lands_initial_gap_above_bottom_apex_symmetric() -> None:
    """Sanity check on the symmetric case: top apex z = bottom apex z + gap."""
    half_x, half_z = SHAPES["sphere"](Hp=220.0)
    z_apex_bot = float(half_z.max())
    _, tvz, _ = top_body_vertices(half_x, half_z, P, PLATE, GAP)
    top_apex_z = float(tvz.min())  # mirrored — apex is the lowest point
    assert top_apex_z == pytest.approx(z_apex_bot + GAP, abs=1e-9)


def test_top_apex_lands_initial_gap_above_bottom_apex_asymmetric() -> None:
    """When sphere is on the bottom and cap on top, the top *cap* apex
    must still sit `initial_gap` above the *sphere* apex."""
    bot_x, bot_z = SHAPES["sphere"](Hp=220.0)
    top_x, top_z = SHAPES["prolate"](Hp=220.0)
    z_apex_bot = float(bot_z.max())
    z_apex_top_profile = float(top_z.max())
    assert z_apex_bot != z_apex_top_profile, (
        "test premise broken — sphere and prolate have the same apex height"
    )

    _, tvz, _ = top_body_vertices(top_x, top_z, P, PLATE, GAP,
                                  bottom_z_apex=z_apex_bot)
    top_apex_z = float(tvz.min())
    assert top_apex_z == pytest.approx(z_apex_bot + GAP, abs=1e-9)


def test_top_plate_above_top_pillar_asymmetric() -> None:
    """Plate sits exactly one pillar-height above the top apex."""
    bot_x, bot_z = SHAPES["sphere"](Hp=220.0)
    top_x, top_z = SHAPES["prolate"](Hp=220.0)
    z_apex_bot = float(bot_z.max())
    z_apex_top_profile = float(top_z.max())
    Hp = -float(top_z[0])

    _, tvz, ttags = top_body_vertices(top_x, top_z, P, PLATE, GAP,
                                      bottom_z_apex=z_apex_bot)
    # The two consecutive vertices tagged "top_disp" form the top edge of
    # the top plate; one of them sits at (0, outer).
    inner = Hp + z_apex_bot + z_apex_top_profile + GAP
    outer = inner + PLATE
    # First vertex in `top_body_vertices` is (0, outer).
    assert tvz[0] == pytest.approx(outer, abs=1e-9)


def test_config_top_profile_optional() -> None:
    """Schema accepts a config with no `top_profile` (legacy) and one with."""
    sym = SimulationConfig.model_validate({
        "geometry": {"pitch_um": P, "Hp_um": 220.0, "plate_um": PLATE,
                     "initial_gap_um": GAP},
        "profile": {"kind": "builtin", "name": "sphere"},
        "solver": {"steps": 10, "disp_um": 10.0},
        "output": {"out_dir": "/tmp/x"},
    })
    assert sym.top_profile is None
    assert sym.effective_top_profile is sym.profile

    asym = SimulationConfig.model_validate({
        "geometry": {"pitch_um": P, "Hp_um": 220.0, "plate_um": PLATE,
                     "initial_gap_um": GAP},
        "profile": {"kind": "builtin", "name": "sphere"},
        "top_profile": {"kind": "builtin", "name": "cap"},
        "solver": {"steps": 10, "disp_um": 10.0},
        "output": {"out_dir": "/tmp/x"},
    })
    assert asym.top_profile is not None
    assert asym.top_profile.name == "cap"
    assert asym.effective_top_profile.name == "cap"


def test_example_config_loads() -> None:
    """The bundled asymmetric example parses cleanly."""
    from pathlib import Path
    path = Path(__file__).resolve().parent.parent / "examples" / "06_asymmetric_pair" / "sphere_vs_cap.yaml"
    cfg = SimulationConfig.from_yaml(path)
    assert cfg.profile.name == "sphere"
    assert cfg.top_profile is not None
    assert cfg.top_profile.name == "cap"
