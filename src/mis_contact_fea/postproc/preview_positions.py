#!/usr/bin/env python3
"""Position preview for all shapes — no simulation, just plot the
initial and final geometric configurations side-by-side.

Edit SHAPE_CONFIGS dict at the top of this file (or pass overrides
via CLI args) and re-run to see updated positions:

    python3 src/preview_positions.py
    python3 src/preview_positions.py --initial-gap-sphere -40 --disp-sphere 35
    open preview_positions.png

Use this to iteratively dial in each shape's INITIAL_GAP and DISP
before committing to an expensive sweep.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mis_contact_fea.mesh.slice_2d import bottom_body_vertices, top_body_vertices
from mis_contact_fea.profiles import SHAPES


# ============================================================
# EDIT THESE — per-shape starting gap and total descent
# ============================================================
SHAPE_CONFIGS = {
    "sphere":     {"initial_gap":  -55.0, "disp":  85.0},
    "oblate":     {"initial_gap":  -30.0, "disp":  40.0},
    "prolate":    {"initial_gap":  -80.0, "disp": 120.0},
    "cone":       {"initial_gap":  -75.0, "disp":  30.9},
    "cap":        {"initial_gap":  -65.0, "disp":  40.0},
    "capcone":    {"initial_gap": -110.0, "disp":  45.0},
    "doublecone": {"initial_gap":  -80.0, "disp":  40.0},
}

PITCH = 198.0
PLATE_THICKNESS = 30.0
Hp = 220.0
# ============================================================


def build_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("preview_positions.png"))
    p.add_argument("--pitch", type=float, default=PITCH)
    p.add_argument("--plate", type=float, default=PLATE_THICKNESS)
    p.add_argument("--Hp", type=float, default=Hp)
    p.add_argument("--shapes", nargs="+", default=list(SHAPE_CONFIGS),
                   help="Subset of shapes to plot.")
    p.add_argument("--dpi", type=int, default=300,
                   help="Output image DPI (higher = sharper, bigger file).")
    p.add_argument("--width-in", type=float, default=14.0,
                   help="Figure width in inches.")
    p.add_argument("--row-in", type=float, default=3.6,
                   help="Per-shape row height in inches.")
    for shape in SHAPE_CONFIGS:
        p.add_argument(f"--initial-gap-{shape}", type=float, default=None,
                       dest=f"initial_gap_{shape}")
        p.add_argument(f"--disp-{shape}", type=float, default=None,
                       dest=f"disp_{shape}")
    return p.parse_args()


def shape_outlines(shape: str, initial_gap: float, disp: float,
                   pitch: float, plate: float, Hp: float):
    half_x, half_z = SHAPES[shape](Hp=Hp)
    bvx, bvz, _ = bottom_body_vertices(half_x, half_z, pitch, plate)
    tvx, tvz, _ = top_body_vertices(half_x, half_z, pitch, plate, initial_gap)
    return np.array(bvx), np.array(bvz), np.array(tvx), np.array(tvz)


def draw_panel(ax, bvx, bvz, tvx, tvz, pitch, title):
    ax.fill(bvx, bvz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.2, zorder=1)
    ax.fill(tvx, tvz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.2, zorder=1)
    ax.axvline(0, linestyle="--", color="#666", linewidth=1, zorder=0)
    ax.axvline(pitch, linestyle="--", color="#666", linewidth=1, zorder=0)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.15)
    # Pin a generous view range so it doesn't auto-fit per panel
    ax.set_xlim(-15, pitch + 15)


def main():
    args = build_args()
    shapes = args.shapes

    # Apply CLI overrides
    for shape in shapes:
        ig_override = getattr(args, f"initial_gap_{shape}", None)
        d_override = getattr(args, f"disp_{shape}", None)
        if ig_override is not None:
            SHAPE_CONFIGS[shape]["initial_gap"] = ig_override
        if d_override is not None:
            SHAPE_CONFIGS[shape]["disp"] = d_override

    n = len(shapes)
    fig, axes = plt.subplots(
        n, 2,
        figsize=(args.width_in, args.row_in * n),
        dpi=args.dpi,
        squeeze=False,
    )

    # Common z-axis range across all panels of one shape so initial vs
    # final are directly comparable.
    for row, shape in enumerate(shapes):
        cfg = SHAPE_CONFIGS[shape]
        ig = cfg["initial_gap"]
        disp = cfg["disp"]

        bvx, bvz, tvx_init, tvz_init = shape_outlines(
            shape, ig, disp, args.pitch, args.plate, args.Hp,
        )
        tvz_final = tvz_init - disp  # top body translated DOWN by DISP

        z_lo = min(bvz.min(), tvz_final.min()) - 10
        z_hi = max(bvz.max(), tvz_init.max()) + 10

        for col, (tvz, lbl) in enumerate([
            (tvz_init,  f"initial   disp = 0"),
            (tvz_final, f"final     disp = {disp:g}"),
        ]):
            ax = axes[row][col]
            draw_panel(
                ax, bvx, bvz, tvx_init, tvz, args.pitch,
                f"{shape} — {lbl}\n"
                f"INITIAL_GAP = {ig:g}   DISP = {disp:g}",
            )
            ax.set_ylim(z_lo, z_hi)
            if col == 0:
                ax.set_ylabel("z (µm)")
            ax.set_xlabel("x (µm)")

    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {args.out}  ({args.dpi} dpi)")
    print()
    print("Current settings:")
    for shape in shapes:
        cfg = SHAPE_CONFIGS[shape]
        print(f"  {shape:<12} INITIAL_GAP={cfg['initial_gap']:>8g}  DISP={cfg['disp']:>8g}")


if __name__ == "__main__":
    main()
