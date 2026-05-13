#!/usr/bin/env python3
"""
Build a two-body contact mesh from STEP parts using Gmsh.

Outputs:
- mesh.xdmf: volume mesh with cell tags (base/structure)
- facet_tags.xdmf: boundary mesh with facet tags (contact + BCs)
- mesh_tags.json: tag ids used by the solver
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import gmsh
import meshio
import numpy as np


def _bbox_for_entities(dim: int, tags: Iterable[int]) -> tuple[float, float, float, float, float, float]:
    bbs = [gmsh.model.getBoundingBox(dim, tag) for tag in tags]
    arr = np.array(bbs, dtype=float)
    mins = arr[:, :3].min(axis=0)
    maxs = arr[:, 3:].max(axis=0)
    return float(mins[0]), float(mins[1]), float(mins[2]), float(maxs[0]), float(maxs[1]), float(maxs[2])


def _select_surfaces_by_z(
    surf_tags: Iterable[int],
    target_z: float,
    tol: float,
    use_max: bool,
) -> list[int]:
    selected: list[int] = []
    for tag in surf_tags:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(2, tag)
        zval = zmax if use_max else zmin
        if abs(zval - target_z) <= tol:
            selected.append(tag)
    return selected


def _select_surfaces_with_fallback(
    surf_tags: Iterable[int],
    target_z: float,
    tol: float,
    use_max: bool,
    z_span: float,
) -> tuple[list[int], float]:
    selected = _select_surfaces_by_z(surf_tags, target_z, tol, use_max)
    if selected:
        return selected, tol
    fallback_tol = max(tol, 0.05 * z_span)
    selected = _select_surfaces_by_z(surf_tags, target_z, fallback_tol, use_max)
    return selected, fallback_tol


def _collect_surface_nodes(surface_tags: Iterable[int]) -> np.ndarray:
    coords_list: list[np.ndarray] = []
    for tag in surface_tags:
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes(2, int(tag))
        if len(node_tags) == 0:
            continue
        coords = np.array(node_coords, dtype=float).reshape(-1, 3)
        coords_list.append(coords)
    if not coords_list:
        return np.zeros((0, 3), dtype=float)
    coords = np.vstack(coords_list)
    return coords


def _min_vertical_gap(base_pts: np.ndarray, struct_pts: np.ndarray) -> float:
    if base_pts.size == 0 or struct_pts.size == 0:
        raise RuntimeError("Cannot compute gap: empty surface point set.")
    try:
        from scipy.spatial import cKDTree  # type: ignore
    except Exception as exc:
        raise RuntimeError("scipy is required for --auto-gap.") from exc

    tree = cKDTree(base_pts[:, :2])
    _, idx = tree.query(struct_pts[:, :2], k=1)
    z_gap = struct_pts[:, 2] - base_pts[idx, 2]
    return float(np.min(z_gap))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a two-body contact mesh from STEP files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-step", type=Path, required=True)
    parser.add_argument("--structure-step", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("mesh_out"))
    parser.add_argument("--mesh-size-mm", type=float, default=2.0)
    parser.add_argument("--clearance-mm", type=float, default=30.0)
    parser.add_argument("--z-offset-mm", type=float, default=0.0)
    parser.add_argument("--x-offset-mm", type=float, default=0.0)
    parser.add_argument("--y-offset-mm", type=float, default=0.0)
    parser.add_argument(
        "--auto-gap",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Adjust z-offset to match the clearance using sampled surface points.",
    )
    parser.add_argument("--surface-tol-mm", type=float, default=0.5)
    parser.add_argument("--mesh-algorithm", type=int, default=10, help="Gmsh 3D algorithm (10=HXT)")
    parser.add_argument("--gmsh-terminal", type=int, default=1)
    parser.add_argument(
        "--flip-xy",
        action="store_true",
        help="Mirror both parts across the XY plane (z -> -z) before alignment.",
    )
    args = parser.parse_args()

    base_step = args.base_step.resolve()
    structure_step = args.structure_step.resolve()
    out_dir = args.out_dir.resolve()

    # Validate input files exist
    if not base_step.exists():
        raise FileNotFoundError(f"Base STEP file not found: {base_step}")
    if not structure_step.exists():
        raise FileNotFoundError(f"Structure STEP file not found: {structure_step}")

    # Validate mesh parameters
    if args.mesh_size_mm <= 0:
        raise ValueError(f"Mesh size must be positive, got {args.mesh_size_mm}")
    if args.clearance_mm < 0:
        raise ValueError(f"Clearance must be non-negative, got {args.clearance_mm}")

    out_dir.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", args.gmsh_terminal)
    gmsh.option.setNumber("Mesh.Algorithm3D", args.mesh_algorithm)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.mesh_size_mm)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.mesh_size_mm)
    gmsh.model.add("contact")

    base_entities = gmsh.model.occ.importShapes(str(base_step))
    struct_entities = gmsh.model.occ.importShapes(str(structure_step))
    gmsh.model.occ.synchronize()

    base_vols = [tag for dim, tag in base_entities if dim == 3]
    struct_vols = [tag for dim, tag in struct_entities if dim == 3]
    if not base_vols or not struct_vols:
        gmsh.finalize()
        raise RuntimeError("STEP import failed: missing volume entities for base/structure.")

    if args.flip_xy:
        gmsh.model.occ.symmetrize(
            [(3, tag) for tag in base_vols + struct_vols],
            0.0,
            0.0,
            1.0,
            0.0,
        )
        gmsh.model.occ.synchronize()
    else:
        gmsh.model.occ.synchronize()
    base_surfaces = sorted(
        {
            abs(tag)
            for dim, tag in gmsh.model.getBoundary([(3, v) for v in base_vols], oriented=False)
            if dim == 2
        }
    )
    struct_surfaces = sorted(
        {
            abs(tag)
            for dim, tag in gmsh.model.getBoundary([(3, v) for v in struct_vols], oriented=False)
            if dim == 2
        }
    )

    base_xmin, base_ymin, base_zmin, base_xmax, base_ymax, base_zmax = _bbox_for_entities(2, base_surfaces)
    base_span = base_zmax - base_zmin
    base_contact, base_contact_tol = _select_surfaces_with_fallback(
        base_surfaces, base_zmax, args.surface_tol_mm, use_max=True, z_span=base_span
    )
    base_fixed, base_fixed_tol = _select_surfaces_with_fallback(
        base_surfaces, base_zmin, args.surface_tol_mm, use_max=False, z_span=base_span
    )
    if not base_contact:
        gmsh.finalize()
        raise RuntimeError(
            "Could not identify base contact surface; increase --surface-tol-mm or check part orientation."
        )
    if not base_fixed:
        gmsh.finalize()
        raise RuntimeError("Could not identify base fixed surface; increase --surface-tol-mm.")
    base_fixed_bbox = _bbox_for_entities(2, base_fixed)
    base_fixed_cx = 0.5 * (base_fixed_bbox[0] + base_fixed_bbox[3])
    base_fixed_cy = 0.5 * (base_fixed_bbox[1] + base_fixed_bbox[4])
    base_fixed_zmin = base_fixed_bbox[2]
    gmsh.model.occ.translate([(3, tag) for tag in base_vols], -base_fixed_cx, -base_fixed_cy, -base_fixed_zmin)
    gmsh.model.occ.synchronize()

    base_surfaces = sorted(
        {
            abs(tag)
            for dim, tag in gmsh.model.getBoundary([(3, v) for v in base_vols], oriented=False)
            if dim == 2
        }
    )
    base_xmin, base_ymin, base_zmin, base_xmax, base_ymax, base_zmax = _bbox_for_entities(2, base_surfaces)
    base_span = base_zmax - base_zmin
    base_contact, base_contact_tol = _select_surfaces_with_fallback(
        base_surfaces, base_zmax, args.surface_tol_mm, use_max=True, z_span=base_span
    )
    base_fixed, base_fixed_tol = _select_surfaces_with_fallback(
        base_surfaces, base_zmin, args.surface_tol_mm, use_max=False, z_span=base_span
    )

    struct_xmin, struct_ymin, struct_zmin, struct_xmax, struct_ymax, struct_zmax = _bbox_for_entities(
        2, struct_surfaces
    )
    struct_span = struct_zmax - struct_zmin
    struct_contact, struct_contact_tol = _select_surfaces_with_fallback(
        struct_surfaces, struct_zmin, args.surface_tol_mm, use_max=False, z_span=struct_span
    )
    struct_top, struct_top_tol = _select_surfaces_with_fallback(
        struct_surfaces, struct_zmax, args.surface_tol_mm, use_max=True, z_span=struct_span
    )

    if not struct_contact:
        gmsh.finalize()
        raise RuntimeError(
            "Could not identify structure contact surface; increase --surface-tol-mm or check part orientation."
        )
    if not struct_top:
        gmsh.finalize()
        raise RuntimeError("Could not identify structure top surface; increase --surface-tol-mm.")

    base_contact_bbox = _bbox_for_entities(2, base_contact)
    base_contact_zmax = base_contact_bbox[5]

    struct_contact_bbox = _bbox_for_entities(2, struct_contact)
    struct_contact_cx = 0.5 * (struct_contact_bbox[0] + struct_contact_bbox[3])
    struct_contact_cy = 0.5 * (struct_contact_bbox[1] + struct_contact_bbox[4])
    struct_contact_zmin = struct_contact_bbox[2]

    dz = (base_contact_zmax - struct_contact_zmin) + float(args.clearance_mm) + float(args.z_offset_mm)
    dx = -struct_contact_cx + float(args.x_offset_mm)
    dy = -struct_contact_cy + float(args.y_offset_mm)
    gmsh.model.occ.translate([(3, tag) for tag in struct_vols], dx, dy, dz)
    gmsh.model.occ.synchronize()

    struct_surfaces = sorted(
        {
            abs(tag)
            for dim, tag in gmsh.model.getBoundary([(3, v) for v in struct_vols], oriented=False)
            if dim == 2
        }
    )
    struct_xmin, struct_ymin, struct_zmin, struct_xmax, struct_ymax, struct_zmax = _bbox_for_entities(
        2, struct_surfaces
    )
    struct_span = struct_zmax - struct_zmin
    struct_contact, struct_contact_tol = _select_surfaces_with_fallback(
        struct_surfaces, struct_zmin, args.surface_tol_mm, use_max=False, z_span=struct_span
    )
    struct_top, struct_top_tol = _select_surfaces_with_fallback(
        struct_surfaces, struct_zmax, args.surface_tol_mm, use_max=True, z_span=struct_span
    )
    if not struct_contact:
        gmsh.finalize()
        raise RuntimeError(
            "Could not identify structure contact surface after translation; increase --surface-tol-mm."
        )
    if not struct_top:
        gmsh.finalize()
        raise RuntimeError(
            "Could not identify structure top surface after translation; increase --surface-tol-mm."
        )

    struct_contact_bbox = _bbox_for_entities(2, struct_contact)
    struct_contact_zmin = struct_contact_bbox[2]
    actual_gap = struct_contact_zmin - base_contact_zmax
    gap_method = "surface_z"

    print(
        "Surface selections:",
        f"base_contact={len(base_contact)}(tol={base_contact_tol:g})",
        f"base_fixed={len(base_fixed)}(tol={base_fixed_tol:g})",
        f"struct_contact={len(struct_contact)}(tol={struct_contact_tol:g})",
        f"struct_top={len(struct_top)}(tol={struct_top_tol:g})",
    )

    if args.auto_gap:
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(2)
        base_pts = _collect_surface_nodes(base_contact)
        struct_pts = _collect_surface_nodes(struct_contact)
        min_gap = _min_vertical_gap(base_pts, struct_pts)
        shift = min_gap - float(args.clearance_mm)
        if abs(shift) > 1e-9:
            gmsh.model.occ.translate([(3, tag) for tag in struct_vols], 0.0, 0.0, -shift)
            gmsh.model.occ.synchronize()
            # After shift, gap is exactly clearance_mm (no need to remesh and recompute)
            actual_gap = float(args.clearance_mm)
        else:
            actual_gap = min_gap
        gap_method = "sampled_xy"
        gmsh.model.mesh.clear()

    tags = {
        "base_volume": 1,
        "structure_volume": 2,
        "base_contact": 11,
        "structure_contact": 12,
        "base_fixed": 13,
        "structure_top": 14,
    }

    gmsh.model.addPhysicalGroup(3, base_vols, tags["base_volume"])
    gmsh.model.setPhysicalName(3, tags["base_volume"], "base")
    gmsh.model.addPhysicalGroup(3, struct_vols, tags["structure_volume"])
    gmsh.model.setPhysicalName(3, tags["structure_volume"], "structure")

    gmsh.model.addPhysicalGroup(2, base_contact, tags["base_contact"])
    gmsh.model.setPhysicalName(2, tags["base_contact"], "base_contact")
    gmsh.model.addPhysicalGroup(2, struct_contact, tags["structure_contact"])
    gmsh.model.setPhysicalName(2, tags["structure_contact"], "structure_contact")
    gmsh.model.addPhysicalGroup(2, base_fixed, tags["base_fixed"])
    gmsh.model.setPhysicalName(2, tags["base_fixed"], "base_fixed")
    gmsh.model.addPhysicalGroup(2, struct_top, tags["structure_top"])
    gmsh.model.setPhysicalName(2, tags["structure_top"], "structure_top")

    # Compute final geometry bounds for pass-through calculation (before finalize)
    final_base_bbox = _bbox_for_entities(2, base_surfaces)
    final_struct_bbox = _bbox_for_entities(2, struct_surfaces)
    base_height = final_base_bbox[5] - final_base_bbox[2]  # zmax - zmin
    struct_height = final_struct_bbox[5] - final_struct_bbox[2]
    # Pass-through displacement = gap + base_height (structure moves past base completely)
    pass_through_mm = float(actual_gap) + base_height

    msh_path = out_dir / "mesh.msh"
    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()

    msh = meshio.read(msh_path)
    if "tetra" not in msh.cells_dict or "triangle" not in msh.cells_dict:
        raise RuntimeError("Expected tetra/triangle cells in mesh output.")

    cell_data = msh.cell_data_dict.get("gmsh:physical", {})
    tetra_tags = cell_data.get("tetra")
    tri_tags = cell_data.get("triangle")
    if tetra_tags is None or tri_tags is None:
        raise RuntimeError("Missing gmsh physical tags in mesh; did physical groups export?")

    mesh = meshio.Mesh(
        points=msh.points,
        cells={"tetra": msh.cells_dict["tetra"]},
        cell_data={"cell_marker": [tetra_tags]},
    )
    facet_mesh = meshio.Mesh(
        points=msh.points,
        cells={"triangle": msh.cells_dict["triangle"]},
        cell_data={"facet_marker": [tri_tags]},
    )

    meshio.write(out_dir / "mesh.xdmf", mesh)
    meshio.write(out_dir / "facet_tags.xdmf", facet_mesh)

    with (out_dir / "mesh_tags.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "tags": tags,
                "units": "mm",
                "clearance_mm": float(args.clearance_mm),
                "z_offset_mm": float(args.z_offset_mm),
                "x_offset_mm": float(args.x_offset_mm),
                "y_offset_mm": float(args.y_offset_mm),
                "actual_gap_mm": float(actual_gap),
                "gap_method": gap_method,
                "base_height_mm": float(base_height),
                "struct_height_mm": float(struct_height),
                "pass_through_mm": float(pass_through_mm),
            },
            f,
            indent=2,
        )

    print(f"Wrote mesh to {out_dir}")
    print(f"Contact gap after alignment: {actual_gap:.3f} mm")
    print(f"Base height: {base_height:.3f} mm, Structure height: {struct_height:.3f} mm")
    print(f"Pass-through displacement: {pass_through_mm:.3f} mm")


if __name__ == "__main__":
    main()
