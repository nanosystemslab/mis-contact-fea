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

from mis_contact_fea.mesh.layout import (
    bottom_body_vertices,
    dedup_consecutive as _dedup_consecutive,
    signed_area as _signed_area,
    top_body_vertices,
)
from mis_contact_fea.profiles import SHAPES


# Re-export for backward compatibility — external callers used to do
# `from mis_contact_fea.mesh.slice_2d import bottom_body_vertices`.
__all__ = ["bottom_body_vertices", "top_body_vertices", "build_mesh", "main"]


# ---------------- gmsh assembly -------------------------------------------


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
    # If main() prepared a profile spec, use the dispatch loader;
    # otherwise fall back to the legacy "look up SHAPES[args.shape]" path.
    profile = getattr(args, "_profile", None)
    top_profile = getattr(args, "_top_profile", None)
    if profile is None:
        half_x, half_z = SHAPES[args.shape](Hp=args.Hp_um)
        top_half_x, top_half_z = half_x, half_z
    else:
        from mis_contact_fea.profiles import load_profile
        half_x, half_z = load_profile(profile, Hp=args.Hp_um)
        if top_profile is None:
            top_half_x, top_half_z = half_x, half_z
        else:
            top_half_x, top_half_z = load_profile(top_profile, Hp=args.Hp_um)

    Rxy_p = float(half_x[1])  # second polyline point is always (Rxy_p, -Hp)
    z_apex_bot = float(np.asarray(half_z).max())
    z_apex_top = float(np.asarray(top_half_z).max())
    asymmetric = top_profile is not None

    bvx, bvz, btags = bottom_body_vertices(half_x, half_z, args.pitch_um, args.plate_um)
    tvx, tvz, ttags = top_body_vertices(
        top_half_x, top_half_z, args.pitch_um, args.plate_um, args.initial_gap_um,
        bottom_z_apex=z_apex_bot if asymmetric else None,
    )

    # Auto-scale guard: warn if mesh_size is much larger than the
    # smallest polyline segment, since gmsh struggles to place a big
    # mesh_size triangle adjacent to a tiny boundary edge. Only clamp
    # when the ratio is severe (>4x), to avoid over-meshing cases
    # where gmsh can stretch boundary edges into larger triangles.
    def _min_seg_len(vx, vz):
        vx = np.asarray(vx, dtype=float)
        vz = np.asarray(vz, dtype=float)
        d = np.hypot(np.diff(vx), np.diff(vz))
        d = d[d > 1e-12]
        return float(d.min()) if len(d) else float("inf")

    min_seg = min(_min_seg_len(bvx, bvz), _min_seg_len(tvx, tvz))
    requested_mesh = float(args.mesh_size_um)
    if requested_mesh > 4.0 * min_seg:
        safe_mesh = max(min_seg * 2.0, 1e-3)
        print(f"  [auto-scale] mesh_size_um={requested_mesh:.2f} much larger than "
              f"min boundary segment {min_seg:.2f}; clamping to {safe_mesh:.2f}")
        args.mesh_size_um = safe_mesh

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
        "top_shape": getattr(args, "top_shape", None) or args.shape,
        "asymmetric": bool(asymmetric),
        "pitch_um": float(args.pitch_um),
        "Hp_um": float(args.Hp_um),
        "plate_um": float(args.plate_um),
        "initial_gap_um": float(args.initial_gap_um),
        "Rxy_p_um": Rxy_p,
        # Back-compat: keep z_apex_um set to the *larger* apex so existing
        # solver auto-displacement heuristics still produce a safe value.
        "z_apex_um": max(z_apex_bot, z_apex_top),
        "z_apex_bot_um": z_apex_bot,
        "z_apex_top_um": z_apex_top,
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
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML simulation config (provides defaults for the flags below).")
    parser.add_argument("--shape", choices=sorted(SHAPES.keys()), default=None,
                        help="Shortcut for --profile-kind builtin --profile-name X.")
    parser.add_argument("--profile-kind", choices=["builtin", "csv", "svg", "step"],
                        default=None, help="Where the polyline comes from.")
    parser.add_argument("--profile-name", default=None,
                        help="Name of builtin profile (sphere, cone, …).")
    parser.add_argument("--profile-path", type=Path, default=None,
                        help="Path to CSV / SVG / STEP file for non-builtin profiles.")
    parser.add_argument("--profile-units", choices=["um", "mm", "cm", "m"], default=None,
                        help="Length unit of values in profile-path (default: um).")
    parser.add_argument("--profile-prepend-pillar", action="store_true",
                        help="Stitch the standard pillar to a loaded CSV/SVG profile "
                             "that contains only the lobe outline.")
    parser.add_argument("--top-profile-kind",
                        choices=["builtin", "csv", "svg", "step"], default=None,
                        help="Override the top body's profile kind (asymmetric pair).")
    parser.add_argument("--top-profile-name", default=None,
                        help="Name of builtin profile for the top body.")
    parser.add_argument("--top-profile-path", type=Path, default=None,
                        help="Path to CSV/SVG/STEP file for the top body's profile.")
    parser.add_argument("--top-profile-units",
                        choices=["um", "mm", "cm", "m"], default=None,
                        help="Length unit of values in --top-profile-path.")
    parser.add_argument("--top-profile-prepend-pillar", action="store_true",
                        help="Stitch the standard pillar to a loaded top-body CSV/SVG.")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--pitch-um", type=float, default=None)
    parser.add_argument("--Hp-um", type=float, default=None,
                        help="Pillar height in micrometres.")
    parser.add_argument("--plate-um", type=float, default=None,
                        help="Plate thickness.")
    parser.add_argument("--initial-gap-um", type=float, default=None,
                        help="Initial z gap between top and bottom structure tips.")
    parser.add_argument("--mesh-size-um", type=float, default=None)
    parser.add_argument("--gmsh-terminal", type=int, default=1)
    args = parser.parse_args()

    # If --config is given, fill any flag the user didn't pass from the
    # config; explicit flags always win over config values. If no config,
    # apply the legacy defaults to maintain backward-compatible CLI.
    if args.config is not None:
        from mis_contact_fea.config import SimulationConfig

        cfg = SimulationConfig.from_yaml(args.config)
        # Profile fields
        if args.profile_kind is None:
            args.profile_kind = cfg.profile.kind
        if args.profile_name is None:
            args.profile_name = cfg.profile.name
        if args.profile_path is None:
            args.profile_path = cfg.profile.path
        if args.profile_units is None:
            args.profile_units = cfg.profile.units
        if not args.profile_prepend_pillar:
            args.profile_prepend_pillar = cfg.profile.prepend_pillar
        # Top-profile fields (None if config doesn't set top_profile)
        if cfg.top_profile is not None:
            if args.top_profile_kind is None:
                args.top_profile_kind = cfg.top_profile.kind
            if args.top_profile_name is None:
                args.top_profile_name = cfg.top_profile.name
            if args.top_profile_path is None:
                args.top_profile_path = cfg.top_profile.path
            if args.top_profile_units is None:
                args.top_profile_units = cfg.top_profile.units
            if not args.top_profile_prepend_pillar:
                args.top_profile_prepend_pillar = cfg.top_profile.prepend_pillar
        # Geometry / mesh fields
        if args.shape is None and cfg.profile.kind == "builtin":
            args.shape = cfg.profile.name
        if args.out_dir is None:
            args.out_dir = cfg.mesh_dir
        if args.pitch_um is None:
            args.pitch_um = cfg.geometry.pitch_um
        if args.Hp_um is None:
            args.Hp_um = cfg.geometry.Hp_um
        if args.plate_um is None:
            args.plate_um = cfg.geometry.plate_um
        if args.initial_gap_um is None:
            args.initial_gap_um = cfg.geometry.initial_gap_um
        if args.mesh_size_um is None:
            args.mesh_size_um = cfg.mesh.characteristic_size_um

    # Legacy defaults for the no-config path.
    if args.profile_kind is None:
        # --shape on its own implies kind=builtin
        args.profile_kind = "builtin"
    if args.profile_units is None:
        args.profile_units = "um"
    if args.shape is None and args.profile_kind == "builtin" and args.profile_name is None:
        args.shape = "sphere"
    if args.profile_kind == "builtin" and args.profile_name is None:
        args.profile_name = args.shape
    if args.out_dir is None:
        args.out_dir = Path("mesh_out_2d")
    if args.pitch_um is None:
        args.pitch_um = 170.0
    if args.Hp_um is None:
        args.Hp_um = 220.0
    if args.plate_um is None:
        args.plate_um = 30.0
    if args.initial_gap_um is None:
        args.initial_gap_um = 10.0
    if args.mesh_size_um is None:
        args.mesh_size_um = 4.0

    # Build a ProfileConfig and stash on args for build_mesh().
    from mis_contact_fea.config import ProfileConfig
    args._profile = ProfileConfig(
        kind=args.profile_kind,
        name=args.profile_name,
        path=args.profile_path,
        units=args.profile_units,
        prepend_pillar=args.profile_prepend_pillar,
    )
    # Asymmetric top body: any --top-profile-* flag triggers it.
    if (args.top_profile_kind is not None
            or args.top_profile_name is not None
            or args.top_profile_path is not None):
        args._top_profile = ProfileConfig(
            kind=args.top_profile_kind or "builtin",
            name=args.top_profile_name,
            path=args.top_profile_path,
            units=args.top_profile_units or "um",
            prepend_pillar=args.top_profile_prepend_pillar,
        )
        args.top_shape = args.top_profile_name or "custom"
    else:
        args._top_profile = None
        args.top_shape = None
    # For shape-name-derived things (output paths, gmsh group names),
    # use either the builtin name or a generic 'custom' label.
    if args.shape is None:
        args.shape = args.profile_name or "custom"

    build_mesh(args)


if __name__ == "__main__":
    main()
