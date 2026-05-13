#!/usr/bin/env python3
"""Per-shape force-displacement panels with contact-onset and release
markers, so it's easy to see where each shape's engagement window
actually starts and ends.

Annotations on each panel:
  • dashed vertical line at the first disp where |Fz| crosses
    --contact-threshold (= contact onset)
  • dashed vertical line at the last  disp where |Fz| crosses
    the threshold (= contact release)
  • the peak (|Fz| max) with its disp
  • text labels with each disp value

Use these markers to choose tighter INITIAL_GAP (start closer to
onset) and DISP (stop right after release).
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SHAPE_COLORS = {
    "sphere":     "#1f77b4",
    "oblate":     "#ff7f0e",
    "prolate":    "#2ca02c",
    "cone":       "#d62728",
    "cap":        "#9467bd",
    "capcone":    "#8c564b",
    "doublecone": "#e377c2",
}


def load_force_csv(path: Path):
    d, f = [], []
    with path.open() as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                try:
                    d.append(float(row[0]))
                    f.append(float(row[1]))
                except ValueError:
                    pass
    return d, f


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs_all_shapes"))
    parser.add_argument(
        "--shapes", nargs="+",
        default=["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"],
    )
    parser.add_argument("--out", type=Path, default=Path("individual_forces.png"))
    parser.add_argument("--abs", action="store_true",
                        help="Plot |Fz| instead of signed Fz.")
    parser.add_argument("--contact-threshold", type=float, default=5.0,
                        help="|Fz| threshold for marking contact onset/release.")
    parser.add_argument("--cols", type=int, default=2,
                        help="Subplot grid columns (rows computed from shape count).")
    args = parser.parse_args()

    # Filter to shapes that actually have data
    data = []
    for s in args.shapes:
        p = args.root / s / "force_displacement.csv"
        if not p.exists():
            continue
        d, f = load_force_csv(p)
        if not d:
            continue
        data.append((s, d, f))

    n = len(data)
    cols = args.cols
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6.5 * cols, 3.8 * rows), dpi=120)
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    print(f"{'shape':<12} {'onset disp':>11} {'peak disp':>11} "
          f"{'release disp':>13} {'window':>9} {'peak |Fz|':>11}")
    print("-" * 72)

    for ax, (shape, d, f) in zip(axes, data):
        f_plot = [abs(v) for v in f] if args.abs else f
        color = SHAPE_COLORS.get(shape, "#222")

        ax.plot(d, f_plot, color=color, linewidth=1.4)
        ax.axhline(0, color="#888", linewidth=0.7, linestyle=":")

        # contact onset = first |Fz| >= threshold
        onset_idx = next(
            (i for i, v in enumerate(f) if abs(v) >= args.contact_threshold), None,
        )
        # contact release = last |Fz| >= threshold
        release_idx = None
        if onset_idx is not None:
            for i in range(len(f) - 1, -1, -1):
                if abs(f[i]) >= args.contact_threshold:
                    release_idx = i
                    break
        # peak = max |Fz|
        peak_idx = max(range(len(f)), key=lambda i: abs(f[i])) if f else None

        if onset_idx is not None:
            d_on = d[onset_idx]
            ax.axvline(d_on, color="#2ca02c", linewidth=1.2, linestyle="--", alpha=0.7)
            ax.text(d_on, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else ax.get_ylim()[1] * 0.05,
                    f"  onset\n  {d_on:.1f}",
                    fontsize=8, color="#2ca02c", va="top", ha="left")
        else:
            d_on = None
        if release_idx is not None:
            d_off = d[release_idx]
            ax.axvline(d_off, color="#d62728", linewidth=1.2, linestyle="--", alpha=0.7)
            ax.text(d_off, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else ax.get_ylim()[1] * 0.05,
                    f"release  \n{d_off:.1f}  ",
                    fontsize=8, color="#d62728", va="top", ha="right")
        else:
            d_off = None
        if peak_idx is not None:
            d_pk = d[peak_idx]
            v_pk = f[peak_idx]
            ax.plot(d_pk, abs(v_pk) if args.abs else v_pk, "o",
                    color="#222", markersize=5, zorder=5)
            ax.annotate(f"peak {abs(v_pk):.0f}\n@ {d_pk:.1f}",
                        (d_pk, abs(v_pk) if args.abs else v_pk),
                        textcoords="offset points", xytext=(5, -5),
                        fontsize=8, color="#222")

        window = (d_off - d_on) if (d_on is not None and d_off is not None) else None
        print(f"{shape:<12} {d_on or 0:>11.2f} {d[peak_idx]:>11.2f} "
              f"{d_off or 0:>13.2f} {window or 0:>9.2f} {abs(f[peak_idx]):>11.1f}")

        ax.set_xlabel("displacement (µm)")
        ax.set_ylabel(("|Fz|" if args.abs else "Fz") + " (µN/µm-depth)")
        ax.set_title(shape, fontsize=11, fontweight="bold")
        ax.grid(alpha=0.25)

    # Hide any unused axes
    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
