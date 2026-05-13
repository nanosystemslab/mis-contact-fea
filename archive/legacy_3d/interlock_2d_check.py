#!/usr/bin/env python3
"""Sanity check: render the notebook's right-half polyline directly and the
full-cell unit cell side-by-side so we can confirm the geometry matches the
notebook. Also show a cone shape where the pillar is fully visible by
construction, for comparison."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from interlock_2d_preview import (
    bottom_body_outline,
    gen_ellipse_axisym,
    top_body_outline,
)

N_ARC = 200


def gen_cone_axisym(Rxy, H, Rxy_p, Hp):
    x = np.array([0.0, Rxy_p, Rxy_p, Rxy, 0.0])
    z = np.array([-Hp, -Hp, 0.0, 0.0, H])
    return x, z


def main():
    fig, axes = plt.subplots(1, 3, figsize=(16, 7), dpi=120)

    # Panel 1: notebook's raw right-half polyline for the sphere
    Rxy, Rz, Rxy_p, Hp = 50.0, 50.0, 25.0, 62.5
    hx, hz = gen_ellipse_axisym(Rxy, Rz, Rxy_p, Hp)
    ax = axes[0]
    ax.fill(hx, hz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.4)
    ax.axvline(0, linestyle="--", linewidth=1.2, color="#888")
    ax.axhline(-Hp, linestyle=":", linewidth=1.0, color="#888")
    ax.axhline(0, linestyle=":", linewidth=1.0, color="#888")
    ax.text(52, -Hp, f"z=-Hp={-Hp:g} (pillar bottom)", fontsize=8, va="center")
    ax.text(52, 0, "z=0 (sphere center)", fontsize=8, va="center")
    z_int = -Rz * np.sqrt(1 - (Rxy_p / Rxy) ** 2)
    ax.axhline(z_int, linestyle=":", linewidth=1.0, color="#888")
    ax.text(52, z_int, f"z={z_int:.1f} (pillar↔sphere merge)", fontsize=8, va="center")
    ax.set_title(
        "Notebook's right-half polyline (sphere)\n"
        f"Hp={Hp:g} but only {-z_int - Hp:.1f}+{Hp + z_int:.1f}={-(z_int):g}-{Hp:g}=... ",
        fontsize=9,
    )
    ax.set_title(
        f"Sphere: notebook right-half polyline\n"
        f"Hp={Hp:g} µm pillar; {-(z_int) - Rz:.1f} µm of shaft visible below sphere\n"
        f"(other {Rz - (-z_int):.0f}+ µm of pillar is buried inside sphere)",
        fontsize=9,
    )
    ax.set_xlim(-5, 130)
    ax.set_ylim(-Hp - 8, Rz + 8)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")
    ax.grid(alpha=0.2)

    # Panel 2: full unit cell with sphere
    ax = axes[1]
    P = 130.0
    plate_thickness = 30.0
    bx, bz = bottom_body_outline(hx, hz, P, plate_thickness)
    tx, tz = top_body_outline(hx, hz, P, plate_thickness, 8.0)
    ax.fill(bx, bz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.4)
    ax.fill(tx, tz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.4)
    ax.axvline(0, linestyle="--", linewidth=1.2, color="#888")
    ax.axvline(P, linestyle="--", linewidth=1.2, color="#888")
    ax.set_title(
        f"Sphere unit cell (P={P:g}, gap=8)\n"
        "Visible pillar shaft is short because the\n"
        "sphere centerline sits at z=0, burying most\n"
        "of the pillar (this matches the notebook).",
        fontsize=9,
    )
    ax.set_xlim(-15, P + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.grid(alpha=0.2)

    # Panel 3: full unit cell with cone (pillar fully visible by construction)
    ax = axes[2]
    Rxy_c, H_c, Rxy_pc, Hp_c = 50.0, 50.0, 25.0, 62.5
    hxc, hzc = gen_cone_axisym(Rxy_c, H_c, Rxy_pc, Hp_c)
    bx, bz = bottom_body_outline(hxc, hzc, P, plate_thickness)
    tx, tz = top_body_outline(hxc, hzc, P, plate_thickness, 8.0)
    ax.fill(bx, bz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.4)
    ax.fill(tx, tz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.4)
    ax.axvline(0, linestyle="--", linewidth=1.2, color="#888")
    ax.axvline(P, linestyle="--", linewidth=1.2, color="#888")
    ax.set_title(
        f"Cone unit cell (same Hp={Hp_c:g})\n"
        "Cone base sits at z=0, so the full Hp pillar\n"
        "is visible underneath. Same pillar dimensions\n"
        "as the sphere panel — just no overlap.",
        fontsize=9,
    )
    ax.set_xlim(-15, P + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.grid(alpha=0.2)

    fig.suptitle(
        "Sanity check: pillar Hp = 62.5 µm in all three panels.\n"
        "For the sphere (panels 1–2), most of the pillar is volumetrically inside the sphere body — the merged outline only shows the shaft below z ≈ −43 µm.\n"
        "For the cone (panel 3), the structure base sits at z=0, so the full pillar is exposed. Same Hp, very different visual.",
        fontsize=10, y=1.04,
    )

    out_path = "interlock_2d_check.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
