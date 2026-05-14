"""Load a right-half axisymmetric profile from an SVG `<path>` element.

Use this when you've sketched a custom lobe shape in Inkscape /
Illustrator / Figma and want to feed it straight to the simulator.

Conventions
-----------
- SVG uses screen coordinates: x grows right, y grows DOWN. We flip
  y to match our z-up convention.
- Units: SVG numbers are in pixels by default. Pass `scale_um_per_px`
  to convert (e.g. if your viewBox is sized so 1 px = 1 µm, set
  `scale_um_per_px=1.0`).
- If your SVG has multiple paths, pass `--path-index N` to pick the
  N-th one (default: 0, the first).
- The resulting polyline must start and end on the symmetry axis
  (x ≈ 0) — we shift it so the leftmost point sits on x=0 before
  validating.

Example workflow
----------------
1. Draw your right-half outline in Inkscape, snapping the start and
   end points to the y-axis (the symmetry axis).
2. Save as plain SVG.
3. Either:
   - Pre-set your Inkscape document units to micrometres, OR
   - Note the px-to-µm scale and pass it in.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def load_svg_profile(
    path: Path | str,
    Hp: float = 220.0,
    units: str = "um",
    samples: int = 200,
    path_index: int = 0,
    flip_y: bool = True,
    prepend_pillar: bool = False,
    Rxy_p: float = 25.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Parse the N-th `<path>` in an SVG file into a right-half polyline.

    Parameters
    ----------
    path : Path or str
        Path to the .svg file.
    Hp : float
        Pillar height in µm (only used when prepend_pillar=True).
    units : str
        Length unit of the SVG coordinate values. One of "um", "mm", "cm", "m".
        Set this to match what one user-coordinate-unit in your SVG represents
        (e.g. if you set up your Inkscape document units as millimetres, use
        "mm"). Values are converted to µm before returning.
    samples : int
        How many points to sample along the SVG path. Higher = smoother
        but more polyline vertices.
    path_index : int
        Which `<path>` element to use if the SVG has multiple.
    flip_y : bool
        Negate y-coordinates so SVG-down maps to FEA-up (default True).
    prepend_pillar, Rxy_p : see load_csv_profile.

    Returns
    -------
    (x, z) : np.ndarray, np.ndarray
        Right-half polyline in µm. The polyline is translated so its
        leftmost x sits on 0 (the symmetry axis).
    """
    from . import um_per_unit

    scale = um_per_unit(units)
    try:
        from svgpathtools import svg2paths
    except ImportError as e:
        raise ImportError(
            "svgpathtools is required for SVG profile loading. "
            "Install with: pip install svgpathtools"
        ) from e

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"profile SVG not found: {path}")

    paths, _attrs = svg2paths(str(path))
    if not paths:
        raise ValueError(f"{path}: no <path> elements found")
    if path_index >= len(paths):
        raise IndexError(
            f"{path}: path_index={path_index} out of range (file has {len(paths)} paths)"
        )

    svg_path = paths[path_index]
    # Uniform-arc-length sampling
    ts = np.linspace(0.0, 1.0, samples)
    pts = np.array([svg_path.point(t) for t in ts])  # complex: real=x, imag=y
    x = pts.real.astype(float)
    y = pts.imag.astype(float)
    if flip_y:
        y = -y

    # Translate so leftmost x is on the axis (x=0).
    x -= x.min()

    # Apply unit scale (input -> µm).
    x *= scale
    z = y * scale

    # Translate so the polyline starts at the axis-tip.
    # SVG paths often start/end at arbitrary points; if the first or
    # last vertex isn't on the axis, snap it.
    if x[0] > 1e-3:
        # Find the axis vertex (smallest x) and rotate the loop.
        i0 = int(np.argmin(x))
        x = np.concatenate([x[i0:], x[:i0]])
        z = np.concatenate([z[i0:], z[:i0]])

    # If first and last aren't both on the axis, the SVG is open / not
    # axisymmetric. Truncate to the segment that lies between axis touches.
    on_axis = np.where(x < 1e-3)[0]
    if len(on_axis) >= 2:
        # Use the segment from the first to the last axis-touching index.
        x = x[on_axis[0]:on_axis[-1] + 1]
        z = z[on_axis[0]:on_axis[-1] + 1]
        # Force the endpoints exactly onto x=0.
        x[0] = 0.0
        x[-1] = 0.0

    if prepend_pillar:
        pillar_x = np.array([0.0, Rxy_p, Rxy_p])
        pillar_z = np.array([-Hp, -Hp, z[0]])
        x = np.concatenate([pillar_x, x])
        z = np.concatenate([pillar_z, z])

    if abs(x[0]) > 1e-3 or abs(x[-1]) > 1e-3:
        raise ValueError(
            f"{path}: polyline endpoints not on axis after parsing "
            f"(start x={x[0]:.4f}, end x={x[-1]:.4f}). The SVG path must "
            f"both start and end on the symmetry axis."
        )

    return x, z
