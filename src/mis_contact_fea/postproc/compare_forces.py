#!/usr/bin/env python3
"""Combined force-displacement plot across all shapes.

Reads each shape's force_displacement.csv from runs_all_shapes/<shape>/
and overlays them on a single figure for direct comparison of
insertion force, peak force, snap-through behavior, etc.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Stable per-shape colors (matching the order the user is used to)
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
    disps, forces = [], []
    with path.open() as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                disps.append(float(row[0]))
                forces.append(float(row[1]))
            except ValueError:
                continue
    return disps, forces


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs_all_shapes"))
    parser.add_argument(
        "--shapes", nargs="+",
        default=["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"],
    )
    parser.add_argument("--out", type=Path, default=Path("force_comparison.png"))
    parser.add_argument("--abs", action="store_true",
                        help="Plot |Fz| instead of signed Fz.")
    parser.add_argument(
        "--mode",
        choices=["raw", "normalized", "contact-aligned"],
        default="raw",
        help="raw: x=displacement (µm). normalized: x=fraction of each shape's "
             "own DISP, so all curves span [0,1]. contact-aligned: x=displacement "
             "minus the disp at which |Fz| first crosses --contact-threshold, "
             "so x=0 marks contact onset for every shape.",
    )
    parser.add_argument(
        "--contact-threshold", type=float, default=5.0,
        help="|Fz| threshold (µN/µm-depth) for detecting contact onset "
             "in contact-aligned mode.",
    )
    parser.add_argument(
        "--trim", action="store_true",
        help="Trim each curve to only the active-contact region: from "
             "the first index where |Fz| crosses --contact-threshold up "
             "to the LAST such index. Implies --mode contact-aligned.",
    )
    args = parser.parse_args()
    if args.trim:
        args.mode = "contact-aligned"

    fig, ax = plt.subplots(figsize=(10, 6.5), dpi=120)

    summary_lines = []

    for shape in args.shapes:
        csv_path = args.root / shape / "force_displacement.csv"
        if not csv_path.exists():
            print(f"skip {shape}: no force_displacement.csv")
            continue
        d, f = load_force_csv(csv_path)
        if not d:
            print(f"skip {shape}: empty CSV")
            continue
        if args.abs:
            f_plot = [abs(v) for v in f]
        else:
            f_plot = f
        if args.mode == "normalized":
            d_max = max(d) if max(d) > 0 else 1.0
            x_plot = [v / d_max for v in d]
        elif args.mode == "contact-aligned":
            # Find first index where |Fz| exceeds the threshold.
            onset_idx = next(
                (i for i, v in enumerate(f) if abs(v) >= args.contact_threshold),
                None,
            )
            if onset_idx is None:
                d_onset = d[0]
                print(f"  {shape}: contact never crosses |Fz|>={args.contact_threshold}; "
                      f"plotting un-shifted")
            else:
                d_onset = d[onset_idx]
                print(f"  {shape}: contact onset at disp={d_onset:.2f} µm "
                      f"(step {onset_idx + 1})")
            x_plot = [v - d_onset for v in d]

            # Optional trim: cut to active-contact region [onset, last-active].
            if args.trim and onset_idx is not None:
                release_idx = next(
                    (i for i, v in enumerate(f[::-1])
                     if abs(v) >= args.contact_threshold),
                    None,
                )
                if release_idx is not None:
                    release_idx = len(f) - 1 - release_idx
                    d_release = d[release_idx]
                    print(f"           release at disp={d_release:.2f} µm "
                          f"(step {release_idx + 1})  → trimmed window "
                          f"{d_release - d_onset:.2f} µm wide")
                    x_plot = x_plot[onset_idx : release_idx + 1]
                    f_plot = f_plot[onset_idx : release_idx + 1]
        else:
            x_plot = d

        peak = max(f_plot, key=abs)
        d_at_peak = d[f_plot.index(peak)]
        color = SHAPE_COLORS.get(shape, None)
        label = f"{shape}  peak={peak:+.0f} µN/µm @ disp={d_at_peak:.0f}"
        ax.plot(x_plot, f_plot, label=label, color=color, linewidth=1.6)

        summary_lines.append((shape, peak, d_at_peak, max(d)))

    ax.axhline(0, color="#888", linewidth=0.8, linestyle="--")
    if args.mode == "contact-aligned":
        ax.axvline(0, color="#444", linewidth=0.9, linestyle=":")
    xlabels = {
        "raw":             "displacement (µm)",
        "normalized":      "normalized displacement (0=initial, 1=final)",
        "contact-aligned": f"displacement past contact onset (µm)  [onset = |Fz|≥{args.contact_threshold:g}]",
    }
    ax.set_xlabel(xlabels[args.mode])
    ax.set_ylabel("reaction force " + ("|Fz|" if args.abs else "Fz") + " (µN/µm-depth)")
    ax.set_title(
        f"Insertion force across {len(summary_lines)} shapes"
        + (" — |Fz|" if args.abs else "")
        + (" — normalized" if args.mode == "normalized" else "")
    )
    ax.legend(fontsize=9, loc="best", framealpha=0.92)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.out}")

    if summary_lines:
        print("\nPeak forces:")
        print(f"  {'shape':<12} {'peak Fz':>14} {'@ disp':>10} {'total disp':>12}")
        for shape, peak, d_at_peak, d_max in summary_lines:
            print(f"  {shape:<12} {peak:>+14.1f} {d_at_peak:>10.1f} {d_max:>12.1f}")


if __name__ == "__main__":
    main()
