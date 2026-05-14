"""Load a right-half axisymmetric profile from a CSV file.

The CSV must have at least two columns interpretable as (x, z) in
micrometres. A header row is auto-detected; common header names like
`x_um, z_um`, `x, z`, `r, z` are tolerated.

The polyline contract matches the built-in shape generators
(`mis_contact_fea.profiles.builtins`):

    - first point on the axis at the bottom of the pillar: (0, -Hp)
    - last point on the axis at the lobe apex:           (0, z_apex)
    - all points in between have x >= 0 (right-half)

If the CSV represents only the lobe outline (not the pillar), set
`prepend_pillar=True` to have the loader stitch the standard pillar
geometry (Rxy_p × Hp) to the bottom of the polyline.

CSV examples
------------
Two-column, with header:

    x_um,z_um
    0,-220
    25,-220
    25,0
    50,0
    49.6,3.1
    ...
    0,50

Two-column, no header (auto-detected):

    0,-220
    25,-220
    ...
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


_HEADER_TOKENS = {"x", "z", "r", "y", "x_um", "z_um", "r_um"}


def _looks_like_header(row: list[str]) -> bool:
    """Return True if the first cell isn't a number — heuristic header detection."""
    try:
        float(row[0])
        return False
    except (ValueError, IndexError):
        return True


def load_csv_profile(
    path: Path | str,
    Hp: float = 220.0,
    units: str = "um",
    prepend_pillar: bool = False,
    Rxy_p: float = 25.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Load (x, z) polyline from a CSV file.

    Parameters
    ----------
    path : Path or str
        Path to the CSV file.
    Hp : float
        Pillar height in µm (only used when `prepend_pillar=True`).
    units : str
        Length unit of the values in the CSV. One of "um", "mm", "cm", "m".
        Values are converted to µm before returning.
    prepend_pillar : bool
        If True, prepend the standard pillar geometry (0,-Hp) -> (Rxy_p,-Hp) ->
        (Rxy_p, z_first) to the loaded polyline. Use this when the CSV
        contains only the lobe outline and you want it stitched to the
        canonical pillar.
    Rxy_p : float
        Pillar shaft radius in µm (only used when `prepend_pillar=True`).

    Returns
    -------
    (x, z) : np.ndarray, np.ndarray
        Right-half polyline in µm, matching the contract in `profiles.builtins`.

    Raises
    ------
    ValueError
        If the polyline doesn't start/end on the axis (x≈0) or has
        non-finite values.
    """
    from . import um_per_unit

    scale = um_per_unit(units)
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"profile CSV not found: {path}")

    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows_raw = [r for r in reader if r and not r[0].startswith("#")]

    if not rows_raw:
        raise ValueError(f"{path}: empty CSV")
    if _looks_like_header(rows_raw[0]):
        rows_raw = rows_raw[1:]
    for r in rows_raw:
        if len(r) < 2:
            raise ValueError(f"{path}: row has fewer than 2 columns: {r!r}")
        rows.append((float(r[0]), float(r[1])))

    x = np.array([p[0] for p in rows], dtype=float) * scale
    z = np.array([p[1] for p in rows], dtype=float) * scale

    if not (np.isfinite(x).all() and np.isfinite(z).all()):
        raise ValueError(f"{path}: non-finite values in polyline")
    if (x < -1e-9).any():
        raise ValueError(f"{path}: x must be >= 0 (right-half profile only)")

    if prepend_pillar:
        # User supplied just the lobe; stitch the standard pillar to the
        # bottom. The lobe's first point should be roughly at the pillar
        # shoulder (Rxy_p, z_int) but we don't enforce it.
        pillar_x = np.array([0.0, Rxy_p])
        pillar_z = np.array([-Hp, -Hp])
        # bridge pillar top -> lobe start if they're not already coincident
        if not (np.isclose(x[0], Rxy_p) and z[0] <= 0):
            pillar_x = np.concatenate([pillar_x, [Rxy_p]])
            pillar_z = np.concatenate([pillar_z, [z[0]]])
        x = np.concatenate([pillar_x, x])
        z = np.concatenate([pillar_z, z])

    # Final contract check.
    if abs(x[0]) > 1e-6:
        raise ValueError(
            f"{path}: first point must be on the axis (x≈0); got x={x[0]:.4f}. "
            f"Set prepend_pillar=True if your CSV only contains the lobe."
        )
    if abs(x[-1]) > 1e-6:
        raise ValueError(
            f"{path}: last point must be on the axis (x≈0); got x={x[-1]:.4f}."
        )

    return x, z
