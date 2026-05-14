"""Load a right-half axisymmetric profile by cross-sectioning a STEP file.

Workflow
--------
1. Load the STEP body via gmsh's OpenCASCADE kernel.
2. Build a cutting plane at y = 0 (default — the body should already be
   centered on the z-axis).
3. Compute the boolean intersection of the body with the half-plane
   x >= 0 (using a fat half-space).
4. Mesh the resulting 1D curves.
5. Walk the boundary CCW to produce a clean (x, z) polyline.

This is the "shortcut" path mentioned in the README: if you have a 3D
axisymmetric STEP file, this gives you a 2D profile usable by the
plane-strain solver without committing to a full 3D contact run.

Assumptions
-----------
- The STEP body is axisymmetric around the z-axis.
- The body is oriented z-up (axis of revolution = z).
- Only one body is in the file (the loader uses the first available
  volume).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def load_step_profile(
    path: Path | str,
    Hp: float = 220.0,
    units: str = "um",
    bbox_pad: float = 100.0,
    cleanup: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a right-half axisymmetric profile from a STEP file.

    Parameters
    ----------
    path : Path or str
        Path to the .step / .stp file.
    Hp : float
        Pillar height in µm (used to extend the section downward if the
        STEP body doesn't already include a pillar). Set Hp=0 to skip.
    units : str
        Length unit of the STEP body geometry. One of "um", "mm", "cm", "m".
        STEP files do carry a unit declaration but they're frequently
        wrong / inconsistent across CAD tools, so we expose this
        explicitly. Values are converted to µm before returning.
    bbox_pad : float
        Padding added to the cutting half-plane, in the input unit.
        Use larger values for bodies bigger than ~1000 input-units.
    cleanup : bool
        Whether to call gmsh.finalize() afterward (default True).

    Returns
    -------
    (x, z) : np.ndarray, np.ndarray
        Right-half polyline in µm. Starts at (0, z_min), ends at (0, z_max).
    """
    from . import um_per_unit

    scale = um_per_unit(units)
    try:
        import gmsh
    except ImportError as e:
        raise ImportError("gmsh is required for STEP profile loading.") from e

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"profile STEP not found: {path}")

    gmsh.initialize(["-noenv"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("step_section")

        # Load the body.
        gmsh.model.occ.importShapes(str(path))
        gmsh.model.occ.synchronize()

        bbox = gmsh.model.getBoundingBox(-1, -1)
        x_min, y_min, z_min, x_max, y_max, z_max = bbox
        x_extent = max(abs(x_max - x_min), bbox_pad)

        # Build a slicing tool: a box that covers x>=0 and the
        # full z range, centered at y=0 (the cutting plane).
        tool_x_size = x_extent + bbox_pad
        tool_y_size = 1.0  # thin slab on y=0
        tool_z_size = max(z_max - z_min, bbox_pad) + 2 * bbox_pad
        tool = gmsh.model.occ.addBox(
            0.0, -tool_y_size / 2, z_min - bbox_pad,
            tool_x_size, tool_y_size, tool_z_size,
        )

        # Body volumes (everything that isn't our tool).
        body_volumes = [
            (3, v) for v in gmsh.model.occ.getEntities(3) if v[1] != tool and isinstance(v, tuple) is False
        ]
        # getEntities returns list of (dim, tag) tuples; rebuild correctly.
        body_volumes = [(d, t) for d, t in gmsh.model.occ.getEntities(3) if t != tool]
        if not body_volumes:
            raise RuntimeError(f"{path}: no volumes loaded from STEP")

        # Intersect each body with the slicing slab; this leaves us with
        # 2D faces on y=0 representing the cross-section.
        sliced, _ = gmsh.model.occ.intersect(
            body_volumes, [(3, tool)], removeObject=True, removeTool=True,
        )
        gmsh.model.occ.synchronize()

        # The intersection of a 3D body with a thin slab is a 3D
        # "shell" — what we actually want is the 2D face on y=0.
        # Extract all surface entities with normals along +/- y.
        section_pts: list[tuple[float, float]] = []
        for dim, tag in gmsh.model.occ.getEntities(2):
            cog = gmsh.model.occ.getCenterOfMass(dim, tag)
            if abs(cog[1]) > 0.5:  # not on the y=0 plane
                continue
            # Walk boundary curves and pull their vertex coordinates.
            for _bdim, btag in gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=False):
                pts = gmsh.model.getValue(1, btag, [0.0])  # start of curve
                if len(pts) >= 3:
                    section_pts.append((pts[0], pts[2]))
                pts = gmsh.model.getValue(1, btag, [1.0])  # end of curve
                if len(pts) >= 3:
                    section_pts.append((pts[0], pts[2]))

        if not section_pts:
            raise RuntimeError(
                f"{path}: STEP section produced no points. Is the body "
                f"actually centered on the z-axis?"
            )

        # Deduplicate + sort by z (the polyline should be roughly monotonic in z).
        pts = np.array(section_pts)
        # Round to dedupe near-duplicates
        rounded = np.round(pts, decimals=3)
        _, idx = np.unique(rounded, axis=0, return_index=True)
        pts = pts[sorted(idx)]
        # Sort by z so we have a clean axis-walk
        order = np.argsort(pts[:, 1])
        pts = pts[order]

        x, z = pts[:, 0] * scale, pts[:, 1] * scale

        # Optional pillar extension (Hp is in µm)
        if Hp > 0 and z.min() > -Hp + 1.0:
            # Find pillar radius at the bottom of what we have
            rxy_p = float(x[np.argmin(z)])
            pillar_x = np.array([0.0, rxy_p])
            pillar_z = np.array([-Hp, -Hp])
            x = np.concatenate([pillar_x, x])
            z = np.concatenate([pillar_z, z])

        # Ensure endpoints touch the axis
        if abs(x[0]) > 1e-3:
            x = np.concatenate([[0.0], x])
            z = np.concatenate([[z[0]], z])
        if abs(x[-1]) > 1e-3:
            x = np.concatenate([x, [0.0]])
            z = np.concatenate([z, [z[-1]]])

        return x, z
    finally:
        if cleanup:
            gmsh.finalize()
