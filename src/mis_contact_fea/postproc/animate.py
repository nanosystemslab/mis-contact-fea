#!/usr/bin/env python3
"""PyVista-based animation of the 2D contact results.

Renders each time step as a cleaned-up frame (proper VTK rendering, no
matplotlib PolyCollection stretching) coloured by von Mises stress on
the displacement-warped mesh.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import csv

import h5py
import meshio
import numpy as np
import pyvista as pv


def _collect_steps(h5_path: Path, function_name: str):
    found: list[tuple[float, str]] = []

    def visitor(name, obj):
        if f"Function/{function_name}/" in name and isinstance(obj, h5py.Dataset):
            tail = name.rsplit("/", 1)[1]
            try:
                step = float(tail.replace("_", "."))
            except ValueError:
                step = len(found)
            found.append((step, name))

    with h5py.File(h5_path, "r") as h5:
        h5.visititems(visitor)
    found.sort(key=lambda kv: kv[0])
    return found


def _find_mesh_paths(h5_path: Path):
    """Locate /Mesh/<name>/geometry and topology paths in a dolfinx XDMF h5."""
    geometry = topology = None

    def visitor(name, obj):
        nonlocal geometry, topology
        if not isinstance(obj, h5py.Dataset):
            return
        low = name.lower()
        if "/mesh/" in low and low.endswith("geometry"):
            geometry = name
        elif "/mesh/" in low and low.endswith("topology"):
            topology = name

    with h5py.File(h5_path, "r") as h5:
        h5.visititems(visitor)
    if geometry is None or topology is None:
        with h5py.File(h5_path, "r") as h5:
            keys = []
            h5.visit(keys.append)
        raise RuntimeError(
            "Couldn't find dolfinx mesh datasets in h5. Contents:\n  "
            + "\n  ".join(keys[:60])
        )
    return geometry, topology


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-dir", type=Path, default=Path("mesh_out_2d"))
    parser.add_argument("--results-dir", type=Path, default=Path("results_2d_test"))
    parser.add_argument("--xdmf-stem", default="contact_results")
    parser.add_argument("--out", type=Path, default=Path("contact_animation_pv.gif"))
    parser.add_argument("--fps", type=float, default=6.0)
    parser.add_argument("--warp-scale", type=float, default=1.0)
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--clim", type=float, nargs=2, default=None,
                        help="Min/max for stress colormap; if unset, auto.")
    parser.add_argument("--clip-percentile", type=float, default=99.0,
                        help="Cap colormap at this percentile to avoid hot-spot saturation.")
    parser.add_argument("--window", type=int, nargs=2, default=[1100, 900])
    parser.add_argument("--show-edges", action="store_true")
    parser.add_argument("--save-frames", action="store_true",
                        help="Also write PNG per frame alongside the GIF.")
    parser.add_argument("--max-frames", type=int, default=60,
                        help="If the simulation has more steps than this, "
                             "subsample evenly to keep the GIF size manageable. "
                             "Set to 0 to use every step.")
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    h5_path = results_dir / f"{args.xdmf_stem}.h5"
    if not h5_path.exists():
        raise FileNotFoundError(h5_path)

    # Optional displacement-vs-step table from force_displacement.csv —
    # row i (after header) gives the total displacement at simulation
    # step i+1.
    disp_per_step = {}
    csv_path = results_dir / "force_displacement.csv"
    if csv_path.exists():
        with csv_path.open() as fh:
            reader = csv.reader(fh)
            next(reader, None)
            for i, row in enumerate(reader, start=1):
                if row:
                    try:
                        disp_per_step[i] = float(row[0])
                    except ValueError:
                        pass

    # Read mesh from the dolfinx h5 itself so the node ordering matches the
    # u / von_mises data. (Reading the gmsh .msh gives a DIFFERENT node
    # order — dolfinx reorders for parallel partitioning — and using that
    # produces totally garbled warped meshes.)
    geom_path, topo_path = _find_mesh_paths(h5_path)
    with h5py.File(h5_path, "r") as h5:
        pts_dx = np.asarray(h5[geom_path])  # (N, gdim)
        tris = np.asarray(h5[topo_path]).astype(np.int64)
    pts3 = pts_dx.astype(np.float64).copy()
    if pts3.shape[1] == 2:
        pts3 = np.column_stack([pts3, np.zeros(len(pts3))])
    elif pts3.shape[1] >= 3:
        pts3 = pts3[:, :3]
    pts3[:, 2] = 0.0
    print(f"mesh: {len(pts3)} nodes, {len(tris)} triangles")

    # Build pyvista UnstructuredGrid: cells stored as flat [n_pts, idx0, idx1, idx2, ...]
    cells = np.column_stack([
        np.full(len(tris), 3, dtype=np.int64),
        tris.astype(np.int64),
    ]).ravel()
    cell_types = np.full(len(tris), pv.CellType.TRIANGLE, dtype=np.uint8)
    grid_ref = pv.UnstructuredGrid(cells, cell_types, pts3)

    u_steps_all = _collect_steps(h5_path, "u")
    vm_steps_all = _collect_steps(h5_path, "von_mises")
    if not u_steps_all or not vm_steps_all:
        raise RuntimeError(f"No u/von_mises datasets in {h5_path}")
    n_total = min(len(u_steps_all), len(vm_steps_all))

    # Subsample evenly for rendering, but keep the FULL sets around for
    # bounds + colormap pre-passes so the camera and color range cover
    # extremes that might fall between subsample points.
    if args.max_frames > 0 and n_total > args.max_frames:
        indices = sorted(set(
            int(round(i)) for i in np.linspace(0, n_total - 1, args.max_frames)
        ))
        u_steps = [u_steps_all[i] for i in indices]
        vm_steps = [vm_steps_all[i] for i in indices]
        n_steps = len(indices)
        print(f"{n_total} total steps; rendering {n_steps} frames "
              f"(every ~{n_total / n_steps:.1f}th step)")
    else:
        u_steps = u_steps_all
        vm_steps = vm_steps_all
        n_steps = n_total
        print(f"{n_steps} time steps (no subsampling)")

    # Pre-pass for colormap range — sample across ALL steps (not just
    # the subsampled set) so peaks between samples still set the scale.
    with h5py.File(h5_path, "r") as h5:
        all_vm = []
        for _, vk in vm_steps_all[:n_total]:
            all_vm.append(np.asarray(h5[vk]).ravel())
        all_vm = np.concatenate(all_vm)
    if args.clim is not None:
        vmin, vmax = args.clim
    else:
        vmax = float(np.percentile(all_vm, args.clip_percentile))
        vmin = 0.0
    print(f"colormap: [{vmin:g}, {vmax:g}] (clip pct={args.clip_percentile})")

    # Pre-pass: compute the union bounding box across ALL displaced frames
    # (using the FULL step list, not the subsample) so the camera is sized
    # to fit motion that might happen between subsample points.
    bx_lo = bx_hi = bz_lo = bz_hi = None
    with h5py.File(h5_path, "r") as h5:
        for i in range(n_total):
            u_arr = np.asarray(h5[u_steps_all[i][1]])
            if u_arr.ndim == 1:
                u_arr = u_arr.reshape(-1, 2)
            elif u_arr.shape[-1] > 2:
                u_arr = u_arr[:, :2]
            x_disp = pts3[:, 0] + args.warp_scale * u_arr[:, 0]
            z_disp = pts3[:, 1] + args.warp_scale * u_arr[:, 1]
            bx_lo = x_disp.min() if bx_lo is None else min(bx_lo, x_disp.min())
            bx_hi = x_disp.max() if bx_hi is None else max(bx_hi, x_disp.max())
            bz_lo = z_disp.min() if bz_lo is None else min(bz_lo, z_disp.min())
            bz_hi = z_disp.max() if bz_hi is None else max(bz_hi, z_disp.max())
    pad = 0.05 * max(bx_hi - bx_lo, bz_hi - bz_lo)
    bx_lo -= pad; bx_hi += pad; bz_lo -= pad; bz_hi += pad
    cx = 0.5 * (bx_lo + bx_hi)
    cz = 0.5 * (bz_lo + bz_hi)
    half_h = 0.5 * (bz_hi - bz_lo)
    half_w = 0.5 * (bx_hi - bx_lo)
    aspect = args.window[0] / args.window[1]
    # parallel_scale is half-height in world units; if the data is wider
    # than the window aspect, expand it so width fits too.
    parallel_scale = max(half_h, half_w / aspect)

    plotter = pv.Plotter(off_screen=True, window_size=tuple(args.window))
    plotter.set_background("white")
    plotter.enable_parallel_projection()
    plotter.camera.focal_point = (cx, cz, 0.0)
    plotter.camera.position = (cx, cz, 1000.0)
    plotter.camera.up = (0.0, 1.0, 0.0)
    plotter.camera.parallel_scale = parallel_scale
    plotter.open_gif(str(args.out), fps=int(args.fps))

    frame_dir = None
    if args.save_frames:
        frame_dir = args.out.parent / (args.out.stem + "_frames")
        frame_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as h5:
        for i in range(n_steps):
            t, uk = u_steps[i]
            _, vk = vm_steps[i]
            u_arr = np.asarray(h5[uk])
            if u_arr.ndim == 1:
                u_arr = u_arr.reshape(-1, 2)
            elif u_arr.shape[-1] > 2:
                u_arr = u_arr[:, :2]
            vm_arr = np.asarray(h5[vk]).ravel()

            warped = grid_ref.copy()
            warped.points = pts3.copy()
            warped.points[:, 0] += args.warp_scale * u_arr[:, 0]
            warped.points[:, 1] += args.warp_scale * u_arr[:, 1]
            warped.cell_data["von_mises"] = vm_arr

            plotter.clear()
            plotter.set_background("white")
            plotter.add_mesh(
                warped,
                scalars="von_mises",
                cmap=args.cmap,
                clim=[vmin, vmax],
                show_edges=args.show_edges,
                edge_color="gray",
                line_width=0.3,
                reset_camera=False,
                scalar_bar_args={
                    "title": "von Mises (MPa)",
                    "color": "black",
                    "title_font_size": 18,
                    "label_font_size": 14,
                    "vertical": True,
                    "n_labels": 5,
                    "fmt": "%.2e",
                    "position_x": 0.88,
                    "position_y": 0.10,
                    "width": 0.05,
                    "height": 0.80,
                },
            )
            disp_label = (
                f"  |  disp = {disp_per_step[int(t)]:.1f} µm"
                if int(t) in disp_per_step else ""
            )
            plotter.add_text(
                f"step {int(t)}/{n_total}  |  frame {i+1}/{n_steps}{disp_label}",
                position="upper_left", color="black", font_size=12,
            )
            # Re-pin the camera each frame because add_mesh / clear can
            # nudge it. enable_parallel_projection + explicit parallel_scale
            # keep the world-to-pixel mapping fixed across frames.
            plotter.camera.focal_point = (cx, cz, 0.0)
            plotter.camera.position = (cx, cz, 1000.0)
            plotter.camera.up = (0.0, 1.0, 0.0)
            plotter.camera.parallel_scale = parallel_scale
            plotter.write_frame()
            if frame_dir is not None:
                plotter.screenshot(frame_dir / f"frame_{i:03d}.png")

    plotter.close()
    print(f"Saved: {args.out}")
    if frame_dir is not None:
        print(f"Per-step PNGs: {frame_dir}/")


if __name__ == "__main__":
    main()
