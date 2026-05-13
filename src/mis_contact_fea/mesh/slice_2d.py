#!/usr/bin/env python3
"""Build the 2D plane-strain unit-cell mesh for the interlock simulation.

Geometry (single unit cell of an offset-by-P/2 Cartesian array):
  - Bottom body: substrate plate spanning x in [0, P], plus pillar +
    structure centered at x=0=P (split across the periodic seam, so
    the right half sits at x in [0, Rxy] and the left half wraps to
    x in [P-Rxy, P]).
  - Top body: substrate plate spanning x in [0, P], plus pillar +
    structure centered at x=P/2 (whole within the cell), mirrored
    vertically and translated so its tip sits initial_gap above the
    bottom body's tip.

Physical groups (1D facets):
  bottom_fixed   - bottom edge of bottom plate (Dirichlet u=0)
  top_disp       - top edge of top plate (Dirichlet u_z = -d, applied)
  bottom_contact - exterior of bottom pillar/structure + exposed plate top
  top_contact    - exterior of top pillar/structure + exposed plate bottom
  sym_left       - x=0 edges of both bodies (Dirichlet u_x=0; equivalent
                   to periodic BC for this mirror-symmetric cell)
  sym_right      - x=P edges of both bodies (Dirichlet u_x=0)

Physical groups (2D cells):
  bottom_volume  - bottom body interior
  top_volume     - top body interior
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gmsh
import meshio
import numpy as np

from mis_contact_fea.profiles import SHAPES


# ---------------- outline construction ------------------------------------

def bottom_body_vertices(half_x, half_z, P, plate_thickness):
    """CCW vertices of the bottom body's outline as parallel arrays.

    Also returns a list `edge_tags` of length len(verts), where
    edge_tags[i] is the physical-group name for the edge from vert i
    to vert (i+1) % len(verts).
    """
    Hp = -float(half_z[0])
    plate_bot = -Hp - plate_thickness
    plate_top = -Hp

    # half_x has N+1 entries (indices 0..N). Original polyline starts at
    # (0, plate_top), ends at (0, z_apex).
    # Mirrored polyline = (P - half_x, half_z), same N+1 entries.
    n_poly = len(half_x)  # = N + 1

    vx, vz, tags = [], [], []

    def add(x, z, edge_tag):
        vx.append(float(x))
        vz.append(float(z))
        tags.append(edge_tag)

    # Plate corners + cell-right-edge cut
    add(0.0,        plate_bot, "bottom_fixed")  # V0
    add(P,          plate_bot, "sym_right")     # V1
    add(P,          plate_top, "sym_right")     # V2
    add(P,          float(half_z[-1]), "bottom_contact")  # V3 = (P, z_apex)

    # Mirrored polyline (curve of pillar half-B's left side, traversed
    # top->bottom). Reversed-then-skip-endpoints: indices N-1 .. 1.
    # This produces vertices V4..V_{N+2}, ending at (P-Rxy_p, plate_top).
    mir_x = (P - np.asarray(half_x))[::-1]
    mir_z = np.asarray(half_z)[::-1]
    for i in range(1, n_poly - 1):
        add(mir_x[i], mir_z[i], "bottom_contact")

    # Implicit edge from V_{N+2}=(P-Rxy_p, plate_top) to next vertex
    # (Rxy_p, plate_top) crosses the exposed plate top; tagged
    # "bottom_contact" by V_{N+2}'s outgoing tag (set in the for-loop above).
    # Now traverse the original polyline indices 1..N-1 (pillar half-A
    # right side, bottom->top, excluding the axis-base at idx 0 and the
    # apex at idx N which we add separately as "sym_left").
    for i in range(1, n_poly - 1):
        add(float(half_x[i]), float(half_z[i]), "bottom_contact")
    # Last polyline point = apex on axis, also where sym_left starts.
    add(float(half_x[-1]), float(half_z[-1]), "sym_left")  # V_{2N+2} = (0, z_apex)

    # Intermediate vertex at (0, plate_top) to break the long sym_left
    # edge into two segments (matches the right side which has (P, plate_top)).
    add(0.0, plate_top, "sym_left")  # V_{2N+3}

    # Closing edge V_{2N+3} -> V0 along x=0 from plate_top to plate_bot
    # is implicit; that edge is tagged "sym_left" by V_{2N+3}'s tag.
    return np.array(vx), np.array(vz), tags


def top_body_vertices(half_x, half_z, P, plate_thickness, initial_gap):
    """CCW vertices of the top body's outline + edge tags."""
    Hp = -float(half_z[0])
    z_apex = float(np.asarray(half_z).max())
    inner = Hp + 2.0 * z_apex + initial_gap  # plate-bottom z (plate sits above)
    outer = inner + plate_thickness

    # Mirrored-in-z polyline of top pillar (in cell coords, centered at x=P/2)
    top_x_right = P / 2.0 + np.asarray(half_x)  # right side of top pillar
    top_z = 2.0 * z_apex + initial_gap - np.asarray(half_z)
    n_poly = len(half_x)

    vx, vz, tags = [], [], []

    def add(x, z, edge_tag):
        vx.append(float(x))
        vz.append(float(z))
        tags.append(edge_tag)

    add(0.0, outer, "top_disp")    # V0
    add(P,   outer, "sym_right")   # V1
    add(P,   inner, "top_contact") # V2 -> (P/2+Rxy_p, inner) below: plate-bottom exposed
    # Right side of top pillar: polyline indices 1..N
    for i in range(1, n_poly):
        add(float(top_x_right[i]), float(top_z[i]), "top_contact")  # V3..V_{N+2}
    # Last point of right side = (P/2, struct_tip).
    # Left side mirrored across x=P/2, reversed, [1:-1]: indices N-1..1
    top_x_left = P / 2.0 - np.asarray(half_x)
    for i in range(n_poly - 2, 0, -1):
        add(float(top_x_left[i]), float(top_z[i]), "top_contact")
    # End of left side at (P/2 - Rxy_p, inner) is the last appended; next
    # we need (0, inner) to traverse the plate-bottom exposed segment.
    add(0.0, inner, "sym_left")
    # Closing edge along x=0 from inner up to outer is sym_left (tagged
    # by the just-appended vertex's outgoing edge).
    return np.array(vx), np.array(vz), tags


# ---------------- gmsh assembly -------------------------------------------

def _signed_area(vx, vz):
    n = len(vx)
    s = 0.0
    for i in range(n):
        j = (i + 1) % n
        s += vx[i] * vz[j] - vx[j] * vz[i]
    return 0.5 * s


def _dedup_consecutive(vx, vz, edge_tags, tol=1e-9):
    """Drop any vertex coincident with its successor (within tol).

    A duplicate creates a zero-length edge in the Gmsh geometry, which
    causes the 2D mesher to hang. The pillar/arc-transition in some of
    the notebook's right-half generators introduces such duplicates.
    """
    n = len(vx)
    keep = []
    for i in range(n):
        j = (i + 1) % n
        if abs(vx[i] - vx[j]) < tol and abs(vz[i] - vz[j]) < tol:
            continue  # drop vertex i; its outgoing edge is zero-length
        keep.append(i)
    return (
        [vx[i] for i in keep],
        [vz[i] for i in keep],
        [edge_tags[i] for i in keep],
    )


def _add_polygon_surface(gmsh, vx, vz, edge_tags, mesh_size):
    """Add a closed polygon to gmsh.geo. Returns (surface_tag, line_tags_by_edge_tag).

    Drops consecutive duplicate vertices (defensive) and auto-orients to
    CCW so Gmsh's plane-surface fill treats the loop's interior as the
    surface domain.
    """
    vx, vz, edge_tags = _dedup_consecutive(list(vx), list(vz), list(edge_tags))

    if _signed_area(vx, vz) < 0:
        vx = list(reversed(vx))
        vz = list(reversed(vz))
        # Each vertex's outgoing tag becomes the incoming tag of the
        # original next vertex; rotate-then-reverse keeps each edge with
        # the same tag.
        edge_tags = list(reversed([edge_tags[(i - 1) % len(edge_tags)] for i in range(len(edge_tags))]))

    point_tags = [
        gmsh.model.geo.addPoint(float(x), float(z), 0.0, mesh_size)
        for x, z in zip(vx, vz)
    ]
    n = len(point_tags)
    line_tags = []
    edges_by_name: dict[str, list[int]] = {}
    for i in range(n):
        j = (i + 1) % n
        line = gmsh.model.geo.addLine(point_tags[i], point_tags[j])
        line_tags.append(line)
        edges_by_name.setdefault(edge_tags[i], []).append(line)
    loop = gmsh.model.geo.addCurveLoop(line_tags)
    surface = gmsh.model.geo.addPlaneSurface([loop])
    return surface, edges_by_name


def build_mesh(args) -> dict:
    shape_fn = SHAPES[args.shape]
    half_x, half_z = shape_fn(Hp=args.Hp_um)

    Rxy_p = float(half_x[1])  # second polyline point is always (Rxy_p, -Hp)

    bvx, bvz, btags = bottom_body_vertices(half_x, half_z, args.pitch_um, args.plate_um)
    tvx, tvz, ttags = top_body_vertices(half_x, half_z, args.pitch_um, args.plate_um, args.initial_gap_um)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", args.gmsh_terminal)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay, robust for 2D
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.mesh_size_um)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.mesh_size_um)
    gmsh.model.add("slice")

    surf_bot, edges_bot = _add_polygon_surface(gmsh, bvx, bvz, btags, args.mesh_size_um)
    surf_top, edges_top = _add_polygon_surface(gmsh, tvx, tvz, ttags, args.mesh_size_um)
    gmsh.model.geo.synchronize()

    tags = {
        "bottom_volume":  1,
        "top_volume":     2,
        "bottom_fixed":  11,
        "top_disp":      12,
        "bottom_contact":13,
        "top_contact":   14,
        "sym_left":      15,
        "sym_right":     16,
    }

    gmsh.model.addPhysicalGroup(2, [surf_bot], tags["bottom_volume"])
    gmsh.model.setPhysicalName(2, tags["bottom_volume"], "bottom_volume")
    gmsh.model.addPhysicalGroup(2, [surf_top], tags["top_volume"])
    gmsh.model.setPhysicalName(2, tags["top_volume"], "top_volume")

    # Helper: union edges-by-name across the two bodies, for tags that
    # appear in both (sym_left, sym_right).
    def union_edges(name):
        return edges_bot.get(name, []) + edges_top.get(name, [])

    for name in ("bottom_fixed", "top_disp", "bottom_contact", "top_contact",
                 "sym_left", "sym_right"):
        lines = union_edges(name)
        if not lines:
            raise RuntimeError(f"No lines assigned to physical group {name!r}")
        gmsh.model.addPhysicalGroup(1, lines, tags[name])
        gmsh.model.setPhysicalName(1, tags[name], name)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    gmsh.model.mesh.generate(2)
    msh_path = out_dir / "mesh.msh"
    gmsh.write(str(msh_path))
    gmsh.finalize()

    # Re-read with meshio to write XDMF + tags JSON (mirrors the 3D pipeline)
    msh = meshio.read(msh_path)
    if "triangle" not in msh.cells_dict or "line" not in msh.cells_dict:
        raise RuntimeError("Expected triangle (2D cells) and line (facets) in mesh.")
    cell_data = msh.cell_data_dict.get("gmsh:physical", {})
    tri_tags = cell_data.get("triangle")
    line_tags = cell_data.get("line")
    if tri_tags is None or line_tags is None:
        raise RuntimeError("Missing gmsh physical tags in mesh.")

    mesh = meshio.Mesh(
        points=msh.points,
        cells={"triangle": msh.cells_dict["triangle"]},
        cell_data={"cell_marker": [tri_tags]},
    )
    facet_mesh = meshio.Mesh(
        points=msh.points,
        cells={"line": msh.cells_dict["line"]},
        cell_data={"facet_marker": [line_tags]},
    )
    meshio.write(out_dir / "mesh.xdmf", mesh)
    meshio.write(out_dir / "facet_tags.xdmf", facet_mesh)

    meta = {
        "tags": tags,
        "units": "um",
        "shape": args.shape,
        "pitch_um": float(args.pitch_um),
        "Hp_um": float(args.Hp_um),
        "plate_um": float(args.plate_um),
        "initial_gap_um": float(args.initial_gap_um),
        "Rxy_p_um": Rxy_p,
        "z_apex_um": float(np.asarray(half_z).max()),
    }
    with (out_dir / "mesh_tags.json").open("w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote mesh to {out_dir}")
    print(f"  triangle cells: {len(tri_tags)}, line facets: {len(line_tags)}")
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the 2D plane-strain unit-cell mesh for the interlock sim.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--shape", choices=sorted(SHAPES.keys()), default="sphere")
    parser.add_argument("--out-dir", type=Path, default=Path("mesh_out_2d"))
    parser.add_argument("--pitch-um", type=float, default=170.0)
    parser.add_argument("--Hp-um", type=float, default=220.0,
                        help="Pillar height in micrometres.")
    parser.add_argument("--plate-um", type=float, default=30.0,
                        help="Plate thickness.")
    parser.add_argument("--initial-gap-um", type=float, default=10.0,
                        help="Initial z gap between top and bottom structure tips.")
    parser.add_argument("--mesh-size-um", type=float, default=4.0)
    parser.add_argument("--gmsh-terminal", type=int, default=1)
    args = parser.parse_args()

    build_mesh(args)


if __name__ == "__main__":
    main()
