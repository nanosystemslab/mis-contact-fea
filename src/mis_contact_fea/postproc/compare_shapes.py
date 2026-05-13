#!/usr/bin/env python3
"""Compose a grid showing initial vs final frame for each shape's sim.

Reads frame PNGs out of runs_all_shapes/<shape>/contact_animation_pv_frames/,
plus the per-shape force CSV to add a tiny "Fz peak" annotation.

Output: comparison_grid.png with rows = shapes, columns = [initial, final].
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


def _peak_fz(csv_path: Path):
    if not csv_path.exists():
        return None
    try:
        lines = csv_path.read_text().strip().splitlines()
        vals = []
        for ln in lines[1:]:
            parts = ln.split(",")
            if len(parts) >= 2:
                vals.append(float(parts[1]))
        if not vals:
            return None
        return max(vals, key=abs)
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs_all_shapes"))
    parser.add_argument(
        "--shapes",
        nargs="+",
        default=["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"],
    )
    parser.add_argument("--out", type=Path, default=Path("comparison_grid.png"))
    args = parser.parse_args()

    rows = []
    for shape in args.shapes:
        frame_dir = args.root / shape / "contact_animation_pv_frames"
        if not frame_dir.exists():
            rows.append((shape, None, None, None))
            continue
        frames = sorted(frame_dir.glob("frame_*.png"))
        if not frames:
            rows.append((shape, None, None, None))
            continue
        first = frames[0]
        last = frames[-1]
        peak = _peak_fz(args.root / shape / "force_displacement.csv")
        rows.append((shape, first, last, peak))

    n = len(rows)
    fig, axes = plt.subplots(n, 2, figsize=(10, 3.2 * n), dpi=120)
    if n == 1:
        axes = axes[None, :]

    for r, (shape, first, last, peak) in enumerate(rows):
        for c, (lbl, p) in enumerate([("initial", first), ("final", last)]):
            ax = axes[r, c]
            if p is None or not p.exists():
                ax.text(0.5, 0.5, f"{shape}\n(no frames)",
                        ha="center", va="center", fontsize=14, color="red",
                        transform=ax.transAxes)
            else:
                img = Image.open(p)
                ax.imshow(img)
            ax.set_xticks([]); ax.set_yticks([])
            if c == 0:
                title = f"{shape}  —  {lbl}"
                if peak is not None:
                    title += f"   |Fz|peak = {abs(peak):.2g} µN/µm"
                ax.set_title(title, fontsize=11, loc="left")
            else:
                ax.set_title(lbl, fontsize=11, loc="left")

    fig.tight_layout()
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
