#!/usr/bin/env python3
"""Scale 2D plane-strain reaction forces to a real-world patch.

Inputs are each shape's force_displacement.csv (Fz in µN per µm of
out-of-plane depth, since the 2D simulation is plane-strain).

For a square lattice with pitch P in both lateral directions (i.e. one
pillar pair per P×P cell), the per-pair force, engagement pressure, and
total force on a finite patch are:

    F_per_pair    = Fz_2D * P                            (µN)
    pressure      = F_per_pair / (P*P) = Fz_2D / P       (µN/µm² = MPa)
    F_total       = pressure * A_patch                   (µN)

The out-of-plane pitch cancels in pressure, so the result depends only
on the in-plane pitch you simulated.
"""
from __future__ import annotations

import argparse
import csv
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
    disps, forces = [], []
    with path.open() as f:
        reader = csv.reader(f)
        next(reader, None)
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
        "--shapes",
        nargs="+",
        default=["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"],
    )
    parser.add_argument("--pitch-um", type=float, default=198.0,
                        help="In-plane lattice pitch (= out-of-plane for square arrays).")
    parser.add_argument("--patch-mm", type=float, nargs=2, default=[10.0, 10.0],
                        help="Patch dimensions in mm (Lx, Ly). Default 10mm × 10mm.")
    parser.add_argument("--out", type=Path, default=Path("patch_forces.png"))
    parser.add_argument("--abs", action="store_true",
                        help="Plot |force| instead of signed.")
    parser.add_argument("--mode", choices=["raw", "contact-aligned"], default="raw")
    parser.add_argument("--contact-threshold", type=float, default=5.0,
                        help="|Fz_2D| threshold for contact-onset detection (µN/µm).")
    args = parser.parse_args()

    P = args.pitch_um
    Lx_um, Ly_um = args.patch_mm[0] * 1000.0, args.patch_mm[1] * 1000.0
    A_patch = Lx_um * Ly_um  # µm²
    n_pairs_per_axis_x = Lx_um / P
    n_pairs_per_axis_y = Ly_um / P
    n_pairs = n_pairs_per_axis_x * n_pairs_per_axis_y

    print(f"Lattice pitch       : {P:.1f} µm  (square array)")
    print(f"Patch size          : {args.patch_mm[0]:g} × {args.patch_mm[1]:g} mm "
          f"= {Lx_um:.0f} × {Ly_um:.0f} µm")
    print(f"Pillar pairs        : {n_pairs_per_axis_x:.1f} × {n_pairs_per_axis_y:.1f} "
          f"≈ {n_pairs:.0f} pairs total")
    print()

    fig, ax = plt.subplots(figsize=(11, 6.8), dpi=120)

    rows = []  # (shape, peak_fz_2d, peak_pressure, peak_total_force, d_at_peak, total_disp)

    for shape in args.shapes:
        csv_path = args.root / shape / "force_displacement.csv"
        if not csv_path.exists():
            print(f"skip {shape}: no force_displacement.csv")
            continue
        d, fz_2d = load_force_csv(csv_path)
        if not d:
            print(f"skip {shape}: empty CSV")
            continue

        # Scale to patch-total force (in newtons for readability)
        # fz_2d is in µN/µm-depth → multiply by patch depth Ly_um (µm) gives
        # the per-column force, then multiply by n_pairs_per_axis_x to get
        # total over the in-plane direction... actually simpler:
        #   F_total = (Fz_2D / P) * A_patch  in µN
        # which works for both signed and abs.
        f_patch_uN = [v * A_patch / P for v in fz_2d]
        f_patch_N = [v * 1e-6 for v in f_patch_uN]  # convert µN → N

        if args.abs:
            f_plot = [abs(v) for v in f_patch_N]
        else:
            f_plot = f_patch_N

        # x-axis handling
        if args.mode == "contact-aligned":
            onset_idx = next(
                (i for i, v in enumerate(fz_2d) if abs(v) >= args.contact_threshold),
                None,
            )
            d_onset = d[onset_idx] if onset_idx is not None else d[0]
            x_plot = [v - d_onset for v in d]
        else:
            x_plot = d

        peak_fz_2d = max(fz_2d, key=abs)
        peak_idx = fz_2d.index(peak_fz_2d)
        peak_pressure_MPa = peak_fz_2d / P  # µN/µm² = MPa (note: signed)
        peak_total_N = peak_fz_2d * A_patch / P * 1e-6
        d_at_peak = d[peak_idx]

        rows.append((
            shape, peak_fz_2d, peak_pressure_MPa, peak_total_N, d_at_peak, max(d),
        ))

        color = SHAPE_COLORS.get(shape)
        label = (f"{shape}  |peak|={abs(peak_total_N):.2f} N "
                 f"({abs(peak_pressure_MPa):.2f} MPa)")
        ax.plot(x_plot, f_plot, label=label, color=color, linewidth=1.6)

    ax.axhline(0, color="#888", linewidth=0.8, linestyle="--")
    if args.mode == "contact-aligned":
        ax.axvline(0, color="#444", linewidth=0.9, linestyle=":")
    ax.set_xlabel("displacement past contact onset (µm)" if args.mode == "contact-aligned"
                  else "displacement (µm)")
    ax.set_ylabel(("|F|" if args.abs else "F") + f" on {args.patch_mm[0]:g}×{args.patch_mm[1]:g} mm patch (N)")
    ax.set_title(
        f"Patch force, pitch P={P:g} µm, patch {args.patch_mm[0]:g}×{args.patch_mm[1]:g} mm"
        f"  ({n_pairs:.0f} pillar pairs)"
    )
    ax.legend(fontsize=9, loc="best", framealpha=0.92)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.out}\n")

    # Summary table
    print(f"  {'shape':<12} {'Fz_2D peak':>14} {'pressure':>12} {'F_patch':>14} {'@ disp':>10}")
    print(f"  {'':<12} {'µN/µm-depth':>14} {'MPa':>12} {'N':>14} {'µm':>10}")
    for shape, fz, p, f_n, d, _ in rows:
        print(f"  {shape:<12} {fz:>+14.1f} {p:>+12.3f} {f_n:>+14.2f} {d:>10.1f}")


if __name__ == "__main__":
    main()
