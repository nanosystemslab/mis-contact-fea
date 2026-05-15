#!/usr/bin/env python3
"""Extract a right-half axisymmetric profile from a STEP file.

The STEP file is assumed to be an *extruded 2D sketch* (a thin slab —
the sketch is the cross-section, and we ignore the small extrusion in
the third direction).

Workflow
--------
1. Load the STEP body via gmsh's OpenCASCADE kernel.
2. Find the largest flat face (= the 2D sketch silhouette).
3. Mesh that face's boundary as 1D line segments.
4. Walk the segments via shared mesh-node tags to recover ordered
   closed cycles. The outer silhouette = the cycle with the largest
   |signed area|.
5. Map the cycle to our right-half profile convention. Two cases:

   - `--axis-along x`: the body's axis-of-revolution is the STEP
     file's x-axis. We take the y >= axis_y portion of the outer
     cycle, map our_x = (their_y - axis_y), our_z = their_x.
   - `--axis-along y`: the axis is STEP's y-axis. We take x >=
     axis_x portion, map our_x = (their_x - axis_x), our_z = their_y.

6. Convert units (STEP→µm) and write a CSV with `x_um,z_um` columns.

Why this is needed
------------------
The package's existing `mis_contact_fea.profiles.from_step` loader
assumes axisymmetric STEPs revolved around z, and slices at y=0. It
doesn't handle the extruded-2D-sketch case where the body's axis of
revolution is x or y. This utility produces a CSV file you can then
feed back through the regular `--profile-kind csv` path (which also
works in the Marimo GUI's CSV upload).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import gmsh
except ImportError:
    sys.stderr.write("Error: gmsh is required (pip install gmsh)\n")
    sys.exit(2)


_UNIT_TO_UM = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1_000_000.0}


def _signed_area(poly: np.ndarray) -> float:
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _find_silhouette_face(extrusion_axis: str) -> int:
    """Return the tag of the largest face flat in the extrusion axis."""
    idx = {"x": (0, 3), "y": (1, 4), "z": (2, 5)}[extrusion_axis]
    candidates = []
    for dim, tag in gmsh.model.occ.getEntities(2):
        bb = gmsh.model.getBoundingBox(dim, tag)
        thickness = bb[idx[1]] - bb[idx[0]]
        if thickness > 0.01:
            continue
        area = gmsh.model.occ.getMass(dim, tag)
        candidates.append((tag, area))
    if not candidates:
        raise RuntimeError(f"no faces flat in {extrusion_axis} found")
    candidates.sort(key=lambda c: -c[1])
    return candidates[0][0]


def _find_cycles(edges: dict[int, set[int]], node_xy: dict[int, tuple[float, float]]
                ) -> list[list[tuple[float, float]]]:
    """Decompose the boundary graph into ordered closed cycles."""
    visited: set[int] = set()
    cycles: list[list[tuple[float, float]]] = []
    for start in list(edges.keys()):
        if start in visited:
            continue
        cyc: list[tuple[float, float]] = []
        prev: int | None = None
        cur = start
        while cur not in visited:
            visited.add(cur)
            cyc.append(node_xy[cur])
            nbrs = [n for n in edges[cur] if n != prev]
            if not nbrs:
                break
            prev, cur = cur, nbrs[0]
            if cur == start:
                break
        if len(cyc) >= 3:
            cycles.append(cyc)
    return cycles


def extract_all_outlines(
    step_path: Path,
    units: str = "mm",
    extrusion_axis: str = "z",
    rev_axis: str = "x",
    axis_offset: float = 0.0,
    mesh_size: float = 0.1,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return the closed outlines of *every* face flat in `extrusion_axis`.

    For a revolved body sliced through its rotation axis, the cross-
    section has two disjoint regions (one on each side of the axis);
    each is its own face in the STEP. This function returns one
    polygon per such face so an SVG writer can paint them all.

    Each polygon is returned as `(x_um, z_um)` in our convention, closed
    (first vertex repeated as last). The z origin is shifted so the
    global minimum z across all polygons sits at 0.
    """
    if rev_axis == extrusion_axis:
        raise ValueError("rev_axis must differ from extrusion_axis")
    radial_axis = ({"x", "y", "z"} - {rev_axis, extrusion_axis}).pop()
    scale = _UNIT_TO_UM[units]

    gmsh.initialize(["-noenv"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("step_to_profile")
        gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()

        # Collect every flat-in-extrusion-axis face. The STEP slab usually
        # has the cross-section repeated at z=0 and z=1; we group faces
        # by their extrusion-axis position and keep only the side with
        # the most geometric content (highest summed face area). That
        # filters out construction artifacts and avoids duplicating the
        # outline.
        from collections import defaultdict
        idx = {"x": (0, 3), "y": (1, 4), "z": (2, 5)}[extrusion_axis]
        ax_i = {"x": 0, "y": 1, "z": 2}
        rev_i = ax_i[rev_axis]
        rad_i = ax_i[radial_axis]
        faces_by_plane: dict[float, list[tuple[int, float]]] = defaultdict(list)
        for dim, tag in gmsh.model.occ.getEntities(2):
            bb = gmsh.model.getBoundingBox(dim, tag)
            if (bb[idx[1]] - bb[idx[0]]) > 0.01:
                continue
            plane = round((bb[idx[0]] + bb[idx[1]]) / 2, 3)
            faces_by_plane[plane].append((tag, gmsh.model.occ.getMass(dim, tag)))
        if not faces_by_plane:
            raise RuntimeError(f"no faces flat in {extrusion_axis} found")
        # Pick the plane (z=0 or z=1) with the most total area.
        best_plane = max(faces_by_plane, key=lambda p: sum(a for _t, a in faces_by_plane[p]))
        unique_faces = [t for t, _a in faces_by_plane[best_plane]]

        # Mesh all boundaries once.
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size / 5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.model.mesh.generate(1)

        nt_all, coords_all, _ = gmsh.model.mesh.getNodes()
        coords_all = np.asarray(coords_all).reshape(-1, 3)
        node_xy = {
            int(t): (float(c[rev_i]), float(c[rad_i]))
            for t, c in zip(nt_all, coords_all)
        }

        outlines: list[tuple[np.ndarray, np.ndarray]] = []
        for face_tag in unique_faces:
            boundary = gmsh.model.getBoundary([(2, face_tag)], oriented=False, recursive=False)
            edges: dict[int, set[int]] = {}
            for _bd, ctag in boundary:
                etypes, _ids, enode_tags = gmsh.model.mesh.getElements(1, ctag)
                for _et, _eids, enodes in zip(etypes, _ids, enode_tags):
                    arr = np.asarray(enodes).reshape(-1, 2)
                    for a, b in arr:
                        a, b = int(a), int(b)
                        edges.setdefault(a, set()).add(b)
                        edges.setdefault(b, set()).add(a)
            cycles = _find_cycles(edges, node_xy)
            if not cycles:
                continue
            outer = max(cycles, key=lambda c: abs(_signed_area(np.asarray(c))))
            poly = np.asarray(outer)
            out_x = (poly[:, 1] - axis_offset) * scale
            out_z = poly[:, 0] * scale
            # Close
            out_x = np.concatenate([out_x, out_x[:1]])
            out_z = np.concatenate([out_z, out_z[:1]])
            outlines.append((out_x, out_z))

        # Shift z so global min is at 0.
        z_min = min(z.min() for _x, z in outlines)
        outlines = [(x, z - z_min) for (x, z) in outlines]
        return outlines
    finally:
        gmsh.finalize()


def extract_outline(
    step_path: Path,
    units: str = "mm",
    extrusion_axis: str = "z",
    rev_axis: str = "x",
    axis_offset: float = 0.0,
    mesh_size: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Single-face wrapper around `extract_all_outlines` for back-compat."""
    outlines = extract_all_outlines(
        step_path, units=units, extrusion_axis=extrusion_axis,
        rev_axis=rev_axis, axis_offset=axis_offset, mesh_size=mesh_size,
    )
    # Pick the polygon with the largest |signed area|.
    return max(outlines, key=lambda p: abs(_signed_area(np.column_stack(p))))


def extract_profile(
    step_path: Path,
    units: str = "mm",
    extrusion_axis: str = "z",
    rev_axis: str = "x",
    axis_offset: float | None = None,
    mesh_size: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a right-half (x_um, z_um) polyline from a STEP file.

    Parameters
    ----------
    extrusion_axis
        Which STEP axis the 2D sketch is extruded along ("x", "y", or "z").
    rev_axis
        Which STEP axis acts as our z (axis of revolution / along-cell axis).
        Must be perpendicular to the extrusion axis.
    axis_offset
        Coordinate of the symmetry axis on the *other* in-plane axis. If
        None, defaults to 0. E.g. with extrusion=z, rev_axis=x, an
        axis_offset of 1.0 means y=1.0 is the symmetry axis; we take the
        portion of the cycle where their_y >= 1.0.
    """
    if rev_axis == extrusion_axis:
        raise ValueError("rev_axis must differ from extrusion_axis")
    radial_axis = ({"x", "y", "z"} - {rev_axis, extrusion_axis}).pop()

    scale = _UNIT_TO_UM[units]
    axis_offset = 0.0 if axis_offset is None else axis_offset

    gmsh.initialize(["-noenv"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("step_to_profile")
        gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()

        silhouette = _find_silhouette_face(extrusion_axis)
        boundary = gmsh.model.getBoundary([(2, silhouette)], oriented=False, recursive=False)
        boundary_tags = [t for _d, t in boundary]

        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size / 5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.model.mesh.generate(1)

        nt_all, coords_all, _ = gmsh.model.mesh.getNodes()
        coords_all = np.asarray(coords_all).reshape(-1, 3)

        ax_i = {"x": 0, "y": 1, "z": 2}
        rev_i = ax_i[rev_axis]
        rad_i = ax_i[radial_axis]

        node_xy = {
            int(t): (float(c[rev_i]), float(c[rev_i] * 0 + c[rad_i]))
            for t, c in zip(nt_all, coords_all)
        }
        # Above: x = rev coord (our z), y = radial coord (our x)

        edges: dict[int, set[int]] = {}
        for ctag in boundary_tags:
            elt_types, _ids, elt_node_tags = gmsh.model.mesh.getElements(1, ctag)
            for _et, _eids, enodes in zip(elt_types, _ids, elt_node_tags):
                arr = np.asarray(enodes).reshape(-1, 2)
                for a, b in arr:
                    a, b = int(a), int(b)
                    edges.setdefault(a, set()).add(b)
                    edges.setdefault(b, set()).add(a)

        cycles = _find_cycles(edges, node_xy)
        if not cycles:
            raise RuntimeError("no boundary cycles found")
        outer = max(cycles, key=lambda c: abs(_signed_area(np.asarray(c))))
        poly = np.asarray(outer)

        # Polygon's "y" coord = radial direction; symmetry axis at axis_offset.
        # Take portion of the loop where (radial >= axis_offset).
        radial_rel = poly[:, 1] - axis_offset
        # Sign with tolerance so floating-point noise on the axis edge
        # (values like -2.7e-15 from a nominal subtract-to-zero) doesn't
        # create spurious sign flips that would be detected as crossings.
        # The on-axis-edge case is handled in the elif branch below.
        radial_tol = 1e-3 * max(np.abs(radial_rel).max(), 1.0)
        signs = np.zeros_like(radial_rel, dtype=int)
        signs[radial_rel > radial_tol] = 1
        signs[radial_rel < -radial_tol] = -1
        signs[signs == 0] = 1  # treat boundary points as upper
        crossings = np.where(np.diff(signs) != 0)[0]

        if len(crossings) >= 2:
            # Standard case: cycle crosses the axis twice → upper-half slice.
            i0, i1 = int(crossings[0]), int(crossings[1])
            slab = poly[i0 + 1:i1 + 2]
            if (slab[:, 1] - axis_offset).mean() < 0:
                slab = np.vstack([poly[i1 + 1:], poly[:i0 + 2]])
            slab = slab.copy()
            slab[0, 1] = axis_offset
            slab[-1, 1] = axis_offset
        elif (radial_rel >= -radial_tol).all():
            # Lobe sits entirely above the axis and (typically) touches
            # the axis along ONE long edge. Trim the closed loop to the
            # three sides that aren't on the axis: short-end-A → lobe
            # surface → short-end-B. Endpoints both snap to axis_offset.
            on_axis = radial_rel < radial_tol
            if on_axis.any() and not on_axis.all():
                # Transitions between on-axis and off-axis stretches.
                trans = np.where(np.diff(on_axis.astype(int)) != 0)[0]
                if len(trans) >= 2:
                    # The two corners bracket the off-axis stretch.
                    # Pick the slice whose interior is OFF the axis.
                    i0, i1 = int(trans[0]), int(trans[1])
                    if on_axis[i0]:
                        # i0+1 is off-axis → take poly[i0+1:i1+2]
                        slab = poly[i0 + 1:i1 + 2]
                    else:
                        # i0+1 is on-axis → wrap around
                        slab = np.vstack([poly[i1 + 1:], poly[:i0 + 2]])
                    slab = slab.copy()
                    slab[0, 1] = axis_offset
                    slab[-1, 1] = axis_offset
                else:
                    slab = poly.copy()
            else:
                slab = poly.copy()
        else:
            # Cycle is entirely below the axis — caller probably picked
            # the wrong axis_offset.
            raise RuntimeError(
                f"outer cycle sits at radial range "
                f"[{poly[:, 1].min():.3f}, {poly[:, 1].max():.3f}] but axis_offset="
                f"{axis_offset:.3f}. Pick an axis_offset between min and max, or "
                f"flip the body."
            )

        out_x = (slab[:, 1] - axis_offset) * scale
        out_z = slab[:, 0] * scale
        if out_z[0] > out_z[-1]:
            out_x = out_x[::-1]
            out_z = out_z[::-1]
        out_z -= out_z.min()  # base at z=0

        # Snap exact-zero endpoints to remove float dust.
        if abs(out_x[0]) < 1e-6 * scale:
            out_x[0] = 0.0
        if abs(out_x[-1]) < 1e-6 * scale:
            out_x[-1] = 0.0

        return out_x, out_z
    finally:
        gmsh.finalize()


def write_csv(path: Path, x: np.ndarray, z: np.ndarray, out_units: str = "um") -> None:
    """Write the polyline as CSV. `out_units` controls both the value
    scale and the header (`x_um,z_um`, `x_mm,z_mm`, etc.).
    """
    factor = 1.0 / _UNIT_TO_UM[out_units]  # µm → out unit
    suffix = out_units
    with path.open("w", encoding="utf-8") as f:
        f.write(f"x_{suffix},z_{suffix}\n")
        for xi, zi in zip(x, z):
            f.write(f"{xi * factor:.6f},{zi * factor:.6f}\n")


def write_svg(path: Path,
              polylines: list[tuple[np.ndarray, np.ndarray]],
              closed: bool = False,
              out_units: str = "um") -> None:
    """Emit an SVG with one `<path>` per polyline. Coordinates are scaled
    to `out_units` so the SVG file is in the same unit the user will pick
    in the GUI (avoiding a double-scale when both the file and the GUI
    dropdown convert to µm)."""
    factor = 1.0 / _UNIT_TO_UM[out_units]  # µm → out unit
    scaled = [(np.asarray(x) * factor, np.asarray(z) * factor) for x, z in polylines]
    all_x = np.concatenate([x for x, _z in scaled])
    all_z = np.concatenate([z for _x, z in scaled])
    x_lo, x_hi = float(all_x.min()), float(all_x.max())
    z_lo, z_hi = float(all_z.min()), float(all_z.max())
    w = max(x_hi - x_lo, 1.0)
    h = max(z_hi - z_lo, 1.0)
    pad = 0.05 * max(w, h)
    vb_w = w + 2 * pad
    vb_h = h + 2 * pad
    vb = f"{x_lo - pad:.3f} {-z_hi - pad:.3f} {vb_w:.3f} {vb_h:.3f}"
    aspect = vb_h / vb_w
    if aspect >= 1:
        px_h, px_w = 600, int(round(600 / aspect))
    else:
        px_w, px_h = 600, int(round(600 * aspect))

    fill = "#cce6ff" if closed else "none"
    paths = []
    for x, z in scaled:
        parts = [f"M {x[0]:.3f} {-z[0]:.3f}"]
        parts.extend(f"L {xi:.3f} {-zi:.3f}" for xi, zi in zip(x[1:], z[1:]))
        if closed:
            parts.append("Z")
        d = " ".join(parts)
        paths.append(
            f'  <path d="{d}" stroke="black" fill="{fill}" '
            f'stroke-width="1" vector-effect="non-scaling-stroke"/>'
        )

    path.write_text(
        f'<?xml version="1.0"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" '
        f'width="{px_w}" height="{px_h}">\n'
        + "\n".join(paths)
        + "\n</svg>\n",
        encoding="utf-8",
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert a STEP file (extruded 2D sketch) to a CSV/SVG of its full 2D cross-section.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("step", type=Path, help="Input .step / .stp file.")
    ap.add_argument("--csv", type=Path, help="Output CSV path (defaults to <step>.csv).")
    ap.add_argument("--svg", type=Path, help="Output SVG path (also written if given).")
    ap.add_argument("--units", choices=list(_UNIT_TO_UM), default="mm",
                    help="Length unit of the STEP body. Values are converted to µm.")
    ap.add_argument("--extrusion-axis", choices=["x", "y", "z"], default="z",
                    help="Axis the 2D sketch is extruded along (the thin direction).")
    ap.add_argument("--rev-axis", choices=["x", "y", "z"], default="x",
                    help="Axis along the body's length (becomes our 'z').")
    ap.add_argument("--axis-offset", type=float, default=0.0,
                    help="Coordinate offset on the radial in-plane axis "
                         "(translates the geometry so x=0 lands wherever you want).")
    ap.add_argument("--mesh-size", type=float, default=0.1,
                    help="Boundary mesh size (in the STEP file's units).")
    ap.add_argument("--right-half", action="store_true",
                    help="Emit a right-half axis-snapped polyline instead of "
                         "the full cross-section. Only works for genuinely "
                         "axisymmetric bodies — the default (full outline) is "
                         "the general case.")
    ap.add_argument("--out-units", choices=list(_UNIT_TO_UM), default=None,
                    help="Unit to write the CSV/SVG values in (default: same "
                         "as --units, so the file is in the natural unit of "
                         "the input STEP).")
    args = ap.parse_args()
    if args.out_units is None:
        args.out_units = args.units

    if not args.step.is_file():
        ap.error(f"STEP file not found: {args.step}")

    csv_path = args.csv if args.csv else args.step.with_suffix(".csv")
    mode = "right-half" if args.right_half else "full-outline"
    print(f"Extracting profile from {args.step.name} (mode={mode}, units={args.units}, "
          f"extrusion={args.extrusion_axis}, rev={args.rev_axis}, axis_offset={args.axis_offset})")
    if not args.right_half:
        polylines = extract_all_outlines(
            args.step, units=args.units,
            extrusion_axis=args.extrusion_axis, rev_axis=args.rev_axis,
            axis_offset=args.axis_offset, mesh_size=args.mesh_size,
        )
        # CSV: dump the largest polygon (single outline; clients that need
        # the full geometry should use the SVG).
        x_csv, z_csv = max(polylines, key=lambda p: abs(_signed_area(np.column_stack(p))))
        print(f"  {len(polylines)} disjoint face(s); largest = {len(x_csv)} pts, "
              f"x [{x_csv.min():.1f}, {x_csv.max():.1f}] µm")
    else:
        x_csv, z_csv = extract_profile(
            args.step, units=args.units,
            extrusion_axis=args.extrusion_axis, rev_axis=args.rev_axis,
            axis_offset=args.axis_offset, mesh_size=args.mesh_size,
        )
        polylines = [(x_csv, z_csv)]
        print(f"  polyline: {len(x_csv)} pts, x [{x_csv.min():.1f}, {x_csv.max():.1f}] µm, "
              f"z [{z_csv.min():.1f}, {z_csv.max():.1f}] µm")

    write_csv(csv_path, x_csv, z_csv, out_units=args.out_units)
    print(f"  wrote {csv_path}")
    if args.svg:
        if args.right_half:
            write_svg(args.svg, [(x_csv, z_csv)], closed=False, out_units=args.out_units)
        else:
            write_svg(args.svg, polylines, closed=True, out_units=args.out_units)
        print(f"  wrote {args.svg}")


if __name__ == "__main__":
    main()
