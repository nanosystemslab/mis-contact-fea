#!/usr/bin/env python3
"""Render the 2D plane-strain mesh produced by build_2d_slice_mesh.py.

Colors triangles by cell tag (bottom_volume vs top_volume) and overlays
facets by tag so we can confirm bottom_fixed / top_disp / contact /
sym_left / sym_right ended up on the right edges.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-dir", type=Path, default=Path("mesh_out_2d"))
    parser.add_argument("--out", type=Path, default=Path("mesh_preview.png"))
    args = parser.parse_args()

    mesh_dir = args.mesh_dir.resolve()
    with (mesh_dir / "mesh_tags.json").open("r") as f:
        meta = json.load(f)
    tags = meta["tags"]
    name_of = {v: k for k, v in tags.items()}

    msh = meshio.read(mesh_dir / "mesh.msh")
    pts = msh.points[:, :2]
    tris = msh.cells_dict["triangle"]
    lines = msh.cells_dict["line"]
    tri_tags = msh.cell_data_dict["gmsh:physical"]["triangle"]
    line_tags = msh.cell_data_dict["gmsh:physical"]["line"]

    fig, ax = plt.subplots(figsize=(11, 13), dpi=120)

    cell_colors = {
        tags["bottom_volume"]: "#CFE0EE",
        tags["top_volume"]:    "#F5DCB8",
    }
    polys = pts[tris]
    facecolors = [cell_colors[int(t)] for t in tri_tags]
    pc = PolyCollection(polys, facecolors=facecolors, edgecolors="#777", linewidths=0.2, zorder=1)
    ax.add_collection(pc)

    facet_styles = {
        tags["bottom_fixed"]:   ("#000000", 3.0, "bottom_fixed (u=0)"),
        tags["top_disp"]:       ("#B22222", 3.0, "top_disp (u_z=-d)"),
        tags["bottom_contact"]: ("#1F4F8B", 2.0, "bottom_contact"),
        tags["top_contact"]:    ("#9C4A14", 2.0, "top_contact"),
        tags["sym_left"]:       ("#888888", 1.2, "sym_left (u_x=0)"),
        tags["sym_right"]:      ("#888888", 1.2, "sym_right (u_x=0)"),
    }
    handles = []
    for tag_id, (color, lw, label) in facet_styles.items():
        mask = line_tags == tag_id
        segs = pts[lines[mask]]
        if len(segs) == 0:
            continue
        lc = LineCollection(segs, colors=color, linewidths=lw, zorder=3)
        ax.add_collection(lc)
        handles.append(plt.Line2D([0], [0], color=color, linewidth=lw, label=f"{label} ({mask.sum()})"))

    ax.set_xlim(pts[:, 0].min() - 10, pts[:, 0].max() + 10)
    ax.set_ylim(pts[:, 1].min() - 15, pts[:, 1].max() + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_title(
        f"2D mesh: {meta['shape']}  |  P={meta['pitch_um']:g} µm  Hp={meta['Hp_um']:g} µm  "
        f"plate={meta['plate_um']:g}  gap={meta['initial_gap_um']:g}\n"
        f"{len(tris)} triangles, {len(lines)} facets",
        fontsize=10,
    )
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        framealpha=0.95,
        borderaxespad=0.0,
    )
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
