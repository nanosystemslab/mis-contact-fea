#!/usr/bin/env python3
"""Read dolfinx XDMF/H5 time series and write a GIF of the deformed mesh
coloured by von Mises stress.

The result file pair is produced by run_2d_slice_sim.py:
    results_2d_test/contact_results.xdmf
    results_2d_test/contact_results.h5

dolfinx writes one h5 dataset per write_function call; this script
discovers them by walking the h5 tree.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import PolyCollection
from PIL import Image


def _collect_function_steps(h5_path: Path, function_name: str):
    """Return list of (step_value, dataset_path) for one Function in the h5 file."""
    found: list[tuple[float, str]] = []

    def visitor(name, obj):
        # dolfinx writes like /Function/u/0  or /Function/u/0.000000  depending on version
        if f"Function/{function_name}/" in name and isinstance(obj, h5py.Dataset):
            tail = name.rsplit("/", 1)[1]
            # tail is the step value (often float-encoded)
            try:
                step = float(tail.replace("_", "."))
            except ValueError:
                step = len(found)
            found.append((step, name))

    with h5py.File(h5_path, "r") as h5:
        h5.visititems(visitor)
    found.sort(key=lambda kv: kv[0])
    return found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-dir", type=Path, default=Path("mesh_out_2d"))
    parser.add_argument("--results-dir", type=Path, default=Path("results_2d_test"))
    parser.add_argument("--xdmf-stem", default="contact_results")
    parser.add_argument("--out", type=Path, default=Path("contact_animation.gif"))
    parser.add_argument("--warp-scale", type=float, default=1.0,
                        help="Multiplier on displacement when warping mesh "
                             "(1 = true displacement, >1 to exaggerate).")
    parser.add_argument("--fps", type=float, default=4.0)
    parser.add_argument("--width", type=int, default=900, help="Frame width in pixels.")
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--clim", type=float, nargs=2, default=None,
                        help="Min/max for stress colormap; if unset, auto-scale "
                             "to global max across all frames.")
    parser.add_argument("--save-frames", action="store_true",
                        help="Also save per-step PNG frames next to the GIF.")
    args = parser.parse_args()

    mesh_dir = args.mesh_dir.resolve()
    results_dir = args.results_dir.resolve()
    h5_path = results_dir / f"{args.xdmf_stem}.h5"
    if not h5_path.exists():
        raise FileNotFoundError(f"No h5 file at {h5_path}")

    msh = meshio.read(mesh_dir / "mesh.msh")
    pts = msh.points[:, :2]
    tris = msh.cells_dict["triangle"]

    u_steps = _collect_function_steps(h5_path, "u")
    vm_steps = _collect_function_steps(h5_path, "von_mises")
    if not u_steps or not vm_steps:
        with h5py.File(h5_path, "r") as h5:
            keys = []
            h5.visit(keys.append)
        raise RuntimeError(
            f"No u/von_mises datasets found in {h5_path}. h5 contents:\n  "
            + "\n  ".join(keys[:50])
        )
    n_steps = min(len(u_steps), len(vm_steps))
    print(f"Found {n_steps} time steps in {h5_path}")

    # Pre-pass: compute global stress range
    vm_max_global = 0.0
    with h5py.File(h5_path, "r") as h5:
        for _, vk in vm_steps[:n_steps]:
            vm = np.asarray(h5[vk]).reshape(-1)
            vm_max_global = max(vm_max_global, float(vm.max()))
    if args.clim is not None:
        vmin, vmax = args.clim
    else:
        vmin, vmax = 0.0, max(vm_max_global, 1e-9)
    print(f"Stress colormap range: [{vmin:g}, {vmax:g}] (MPa)")

    # Frame bounds
    x_lo, x_hi = pts[:, 0].min(), pts[:, 0].max()
    z_lo, z_hi = pts[:, 1].min(), pts[:, 1].max()
    margin = 0.05 * max(x_hi - x_lo, z_hi - z_lo)
    aspect = (z_hi - z_lo) / (x_hi - x_lo)
    fig_w_in = args.width / 150
    fig_h_in = fig_w_in * aspect * 1.05 + 0.6

    frames = []
    with h5py.File(h5_path, "r") as h5:
        for i in range(n_steps):
            t, uk = u_steps[i]
            _, vk = vm_steps[i]
            u_arr = np.asarray(h5[uk])
            if u_arr.ndim == 1:
                u_arr = u_arr.reshape(-1, 2)
            elif u_arr.shape[-1] > 2:
                u_arr = u_arr[:, :2]
            vm_arr = np.asarray(h5[vk]).reshape(-1)

            pts_warp = pts + args.warp_scale * u_arr
            polys = pts_warp[tris]

            fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in), dpi=150)
            pc = PolyCollection(
                polys, array=vm_arr,
                cmap=args.cmap, edgecolor="#888", linewidth=0.1,
            )
            pc.set_clim(vmin, vmax)
            ax.add_collection(pc)
            cbar = fig.colorbar(pc, ax=ax, shrink=0.85, label="von Mises (MPa)")
            ax.set_xlim(x_lo - margin, x_hi + margin)
            ax.set_ylim(z_lo - margin, z_hi + margin)
            ax.set_aspect("equal")
            ax.set_xlabel("x (µm)")
            ax.set_ylabel("z (µm)")
            ax.set_title(f"step {i+1}/{n_steps}  (t = {t:g})", fontsize=10)
            fig.tight_layout()

            # Render to memory
            fig.canvas.draw()
            buf = np.asarray(fig.canvas.buffer_rgba())
            frames.append(Image.fromarray(buf).convert("RGB"))
            plt.close(fig)

    duration_ms = int(round(1000.0 / args.fps))
    # disposal=2 ("restore to background") + optimize=False prevents the
    # "smearing" artifact where successive frames overlay each other
    # instead of replacing.
    frames[0].save(
        args.out,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
        disposal=2,
    )
    print(f"Saved: {args.out}  ({len(frames)} frames @ {args.fps:g} fps)")

    if args.save_frames:
        out_dir = args.out.parent / (args.out.stem + "_frames")
        out_dir.mkdir(parents=True, exist_ok=True)
        for i, frame in enumerate(frames):
            frame.save(out_dir / f"frame_{i:03d}.png")
        print(f"Saved {len(frames)} per-step PNGs to {out_dir}/")


if __name__ == "__main__":
    main()
