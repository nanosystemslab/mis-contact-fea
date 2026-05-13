#!/usr/bin/env python3
"""Visualize the 2D plane-strain unit cell for the interlock simulation.

Bottom body: substrate plate + pillar centered on the periodic seam (x=0=P),
so its right half sits at x in [0, Rxy] and its left half wraps to
x in [P-Rxy, P].
Top body: substrate plate (above) + pillar centered at x=P/2, mirrored in z.
Top body is laterally offset by P/2 from bottom body, which is what allows
the lobes to interlock when pushed together.
Lateral edges (x=0, x=P) take Dirichlet u_x=0; for this offset-by-P/2
layout the cell is mirror-symmetric, so this is exactly equivalent to a
true periodic BC.

Each body is rendered as a single closed outline (no overlapping fills),
so internal pillar/plate joins don't show as visible seams.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

N_ARC = 200


def gen_ellipse_axisym(Rxy, Rz, Rxy_p, Hp, n=N_ARC):
    z_int = -Rz * np.sqrt(1 - (Rxy_p / Rxy) ** 2)
    pillar_x = np.array([0.0, Rxy_p, Rxy_p])
    pillar_z = np.array([-Hp, -Hp, z_int])
    theta_start = np.arctan2(z_int / Rz, Rxy_p / Rxy)
    theta_end = np.pi / 2
    t = np.linspace(theta_start, theta_end, n)
    arc_x = Rxy * np.cos(t)
    arc_z = Rz * np.sin(t)
    return np.concatenate([pillar_x, arc_x]), np.concatenate([pillar_z, arc_z])


def bottom_body_outline(half_x, half_z, P, plate_thickness):
    """Single closed outline for the bottom body (CCW traversal):
    plate corners → up cell-right edge through pillar half-B's axis →
    pillar half-B's left profile (down to plate top) → exposed plate top
    → pillar half-A's right profile (up to apex) → down cell-left edge
    through pillar half-A's axis → close along plate-left edge.
    """
    Hp = -half_z[0]
    plate_bot = -Hp - plate_thickness
    plate_top = -Hp
    z_apex = float(half_z.max())

    pts_x = [0.0, P, P, P]
    pts_z = [plate_bot, plate_bot, plate_top, z_apex]

    # Pillar half-B (centered at x=P): mirror across x=P, traverse top→base.
    # In our right-half polyline, index 0 = axis-base, N = axis-apex.
    # Reversed-then-skip-endpoints traverses just-below-apex → just-above-base.
    mirB_x = (P - half_x)[::-1]
    mirB_z = half_z[::-1]
    pts_x.extend(mirB_x[1:-1].tolist())
    pts_z.extend(mirB_z[1:-1].tolist())

    # Endpoint of half-B's left profile is at (P-Rxy_p, plate_top).
    # Implicit segment to half-A's right base at (Rxy_p, plate_top).
    # Pillar half-A (centered at x=0): traverse base → apex on the right side.
    pts_x.extend(half_x[1:].tolist())
    pts_z.extend(half_z[1:].tolist())

    # Down cell-left edge through half-A's axis to plate_top, then plate_bot.
    pts_x.append(0.0)
    pts_z.append(plate_bot)

    return np.array(pts_x), np.array(pts_z)


def top_body_outline(half_x, half_z, P, plate_thickness, initial_gap):
    """Single closed outline for the top body (CCW traversal).

    Top body's pillar is whole within the cell at x_c = P/2; structure tip
    points DOWN. The polyline is mirrored in z (z -> 2*z_apex_bot + gap - z)
    and shifted in x by +P/2.
    """
    Hp = -half_z[0]
    z_apex = float(half_z.max())
    top_plate_inner = Hp + 2 * z_apex + initial_gap
    top_plate_outer = top_plate_inner + plate_thickness

    top_x_right = P / 2 + half_x
    top_z = 2 * z_apex + initial_gap - half_z

    pts_x = [0.0, P, P]
    pts_z = [top_plate_outer, top_plate_outer, top_plate_inner]

    # Right side of pillar: polyline index 1..N (skip axis-base at idx 0).
    pts_x.extend(top_x_right[1:].tolist())
    pts_z.extend(top_z[1:].tolist())

    # Left side: mirror x across P/2, reverse, skip apex (just visited)
    # and axis-base (interior).
    top_x_left = P / 2 - half_x
    pts_x.extend(top_x_left[::-1][1:-1].tolist())
    pts_z.extend(top_z[::-1][1:-1].tolist())

    pts_x.append(0.0)
    pts_z.append(top_plate_inner)

    return np.array(pts_x), np.array(pts_z)


def draw_panel(ax, Rxy, Rz, Rxy_p, Hp, P, plate_thickness, initial_gap, title):
    half_x, half_z = gen_ellipse_axisym(Rxy, Rz, Rxy_p, Hp)
    z_apex = float(half_z.max())

    bx, bz = bottom_body_outline(half_x, half_z, P, plate_thickness)
    tx, tz = top_body_outline(half_x, half_z, P, plate_thickness, initial_gap)

    ax.fill(bx, bz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.4, zorder=1)
    ax.fill(tx, tz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.4, zorder=1)

    ax.axvline(0, linestyle="--", linewidth=1.4, color="#555", zorder=0)
    ax.axvline(P, linestyle="--", linewidth=1.4, color="#555", zorder=0)

    plate_bot_z = -Hp - plate_thickness
    top_outer = Hp + 2 * z_apex + initial_gap + plate_thickness

    ax.annotate(
        "",
        xy=(P / 2, top_outer + 1),
        xytext=(P / 2, top_outer + 14),
        arrowprops=dict(arrowstyle="->", color="#B22222", lw=2),
    )
    ax.text(
        P / 2, top_outer + 16, r"$u_z = -d$",
        ha="center", va="bottom", fontsize=10, color="#B22222",
    )

    for x in np.linspace(2, P - 2, 14):
        ax.plot([x, x - 4], [plate_bot_z, plate_bot_z - 5],
                color="black", linewidth=0.9, zorder=2)
    ax.text(
        P / 2, plate_bot_z - 13, r"$u_x = u_z = 0$ (fixed)",
        ha="center", va="top", fontsize=9,
    )

    mid_z = (top_outer + plate_bot_z) / 2
    ax.text(-12, mid_z, "$u_x = 0$\n(sym ≡ periodic)",
            ha="center", va="center", fontsize=8, color="#333", rotation=90)
    ax.text(P + 12, mid_z, "$u_x = 0$\n(sym ≡ periodic)",
            ha="center", va="center", fontsize=8, color="#333", rotation=-90)

    pitch_y = top_outer + 30
    ax.annotate(
        "", xy=(P, pitch_y), xytext=(0, pitch_y),
        arrowprops=dict(arrowstyle="<->", color="#222", lw=1.2),
    )
    ax.text(P / 2, pitch_y + 2, f"pitch P = {P:g} µm",
            ha="center", va="bottom", fontsize=9)

    ax.set_xlim(-30, P + 30)
    ax.set_ylim(plate_bot_z - 28, pitch_y + 18)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.15)


def main():
    Rxy, Rz, Rxy_p = 50.0, 50.0, 25.0
    Hp = 220.0  # tall enough that top cap can fully pass below bottom cap
    P = 170.0   # 15 µm lobe overlap at peak interference (gentle snap)
    plate_thickness = 30.0

    # Engagement geometry (with offset P/2):
    #   - Top sphere center descends from z = z_top_init toward z = 0 then below.
    #   - At z_top_center = 0 (centers aligned at the equator): "peak
    #     lateral interference", top tip at z = -Rz, lobes overlap by
    #     (2*Rxy - P/2). This is the snap-in resistance peak.
    #   - "Full engagement" (cap entirely past mate): top cap's TOP
    #     reaches the level of bottom cap's BOTTOM, i.e.
    #     z_top_tip + 2*Rz <= -Rz → z_top_tip <= -3*Rz.
    #
    # In this script's parameterization, top tip z = +Rz + initial_gap, so
    # initial_gap = -2*Rz puts top tip at z=-Rz  (peak interference)
    # initial_gap = -4*Rz puts top tip at z=-3*Rz (full engagement: top
    #               cap's top edge touches bottom cap's bottom)
    peak_gap = -2.0 * Rz       # -100; top tip at z=-50
    full_engage_gap = -4.0 * Rz  # -200; top tip at z=-150
    overlap_at_peak = 2.0 * Rxy - P / 2.0
    plate_clear_peak = Hp - Rz                # = 170 with Hp=220
    plate_clear_full_engage = Hp - 3.0 * Rz   # =  70 with Hp=220

    fig, axes = plt.subplots(1, 3, figsize=(15, 12), dpi=120)
    draw_panel(
        axes[0], Rxy, Rz, Rxy_p, Hp, P, plate_thickness,
        initial_gap=10.0,
        title="(a) Initial (10 µm tip gap)",
    )
    draw_panel(
        axes[1], Rxy, Rz, Rxy_p, Hp, P, plate_thickness,
        initial_gap=peak_gap,
        title=(
            f"(b) Peak interference (centers at equator)\n"
            f"lobe overlap {overlap_at_peak:g} µm, plate clear {plate_clear_peak:g} µm"
        ),
    )
    draw_panel(
        axes[2], Rxy, Rz, Rxy_p, Hp, P, plate_thickness,
        initial_gap=full_engage_gap,
        title=(
            f"(c) Full engagement (top cap fully below bottom cap)\n"
            f"top tip at z=-{3*Rz:g}, plate clear {plate_clear_full_engage:g} µm"
        ),
    )

    fig.suptitle(
        f"2D plane-strain unit cell — sphere R={Rz:g} µm, pillar Hp={Hp:g} µm, "
        f"plate={plate_thickness:g} µm, pitch P={P:g} µm (offset P/2)",
        fontsize=12, y=1.02,
    )

    out_path = "interlock_2d_layout.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
