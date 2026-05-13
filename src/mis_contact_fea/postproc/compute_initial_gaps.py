#!/usr/bin/env python3
"""For each lobe shape, find initial_gap such that the closest-point
Euclidean distance between the bottom right-half lobe (centered at
x=0=P, on the periodic seam) and the top lobe (centered at x=P/2,
mirrored vertically) equals a target value.

We baseline the target at "sphere at initial_gap=-50 with the current
pitch", so all shapes start visually as separated as the sphere does.

Usage:
    python3 src/compute_initial_gaps.py --pitch 198 --baseline-shape sphere \
        --baseline-gap -50
    # Or override the target distance directly:
    python3 src/compute_initial_gaps.py --pitch 198 --target-um 6.78
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from mis_contact_fea.profiles import SHAPES


def densify_polyline(x, z, n_per_segment=8):
    """Insert intermediate points along each segment so closest-point sampling
    isn't blind to the middle of long edges."""
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    out_x, out_z = [], []
    for i in range(len(x) - 1):
        ts = np.linspace(0, 1, n_per_segment, endpoint=False)
        out_x.extend(x[i] + (x[i+1] - x[i]) * ts)
        out_z.extend(z[i] + (z[i+1] - z[i]) * ts)
    out_x.append(x[-1])
    out_z.append(z[-1])
    return np.array(out_x), np.array(out_z)


def lobe_points_bot(half_x, half_z, pitch):
    """Return (x, z) sampled points along the bottom body's two lobe halves
    inside the cell, skipping the axis-base segment that's interior to the
    plate-merged body."""
    # half[0] is axis-base at (0, -Hp); drop it. half[1] onward is the surface.
    hx, hz = densify_polyline(half_x[1:], half_z[1:])
    # Right half (centered at x=0): x in [0, Rxy]
    right = np.column_stack([hx, hz])
    # Left half (centered at x=P, wraps to x in [P-Rxy, P])
    left = np.column_stack([pitch - hx, hz])
    return np.vstack([right, left])


def lobe_points_top(half_x, half_z, pitch, initial_gap):
    """Return (x, z) sampled points along the top body's lobe (centered at
    x = pitch/2, mirrored in z)."""
    z_apex = float(np.asarray(half_z).max())
    hx, hz = densify_polyline(half_x[1:], half_z[1:])
    top_z = 2.0 * z_apex + initial_gap - hz
    right = np.column_stack([pitch / 2.0 + hx, top_z])
    left = np.column_stack([pitch / 2.0 - hx, top_z])
    return np.vstack([right, left])


def closest_distance(half_x, half_z, pitch, initial_gap):
    bot = lobe_points_bot(half_x, half_z, pitch)
    top = lobe_points_top(half_x, half_z, pitch, initial_gap)
    # pairwise distances; for ~80*8=640 pts × 640 pts this is ~400k flops, fine.
    dx = top[:, None, 0] - bot[None, :, 0]
    dz = top[:, None, 1] - bot[None, :, 1]
    d = np.sqrt(dx * dx + dz * dz)
    return float(d.min())


def find_gap_for_target(half_x, half_z, pitch, target,
                        gap_max=300.0, gap_min=-300.0, n_samples=300, tol=1e-3,
                        max_iter=80):
    """Return the LARGEST initial_gap such that the closest-point distance
    between the two bodies equals `target`.

    closest_distance(initial_gap) is non-monotonic — it decreases as the top
    body descends, reaches a minimum at peak engagement, then increases again
    after the bodies pass each other. We want the "approaching" solution
    (higher initial_gap, bodies still separated and about to engage), not the
    "departing" mirror solution.
    """
    gaps = np.linspace(gap_max, gap_min, n_samples)
    dists = np.array([closest_distance(half_x, half_z, pitch, float(g)) for g in gaps])

    # First index where the distance has dropped to <= target. Because we
    # sample from high gap (far apart) to low gap (closer or past), this
    # picks the approaching branch.
    below = dists <= target
    if not below.any():
        # Even at peak engagement the bodies never get this close; return
        # the gap that gives the smallest distance achievable.
        return float(gaps[dists.argmin()])
    first = int(below.argmax())
    if first == 0:
        # The very first sample is already below target — bracket failed.
        # Push gap_max higher or return endpoint.
        return float(gaps[0])

    hi = float(gaps[first - 1])  # this one had dist > target
    lo = float(gaps[first])      # this one has dist <= target
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        d = closest_distance(half_x, half_z, pitch, mid)
        if abs(d - target) < tol:
            return mid
        if d > target:
            hi = mid          # still too far, need MORE descent (smaller gap)
        else:
            lo = mid          # got too close, need LESS descent (larger gap)
    return 0.5 * (lo + hi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pitch", type=float, default=198.0)
    parser.add_argument("--baseline-shape", default="sphere")
    parser.add_argument("--baseline-gap", type=float, default=-50.0)
    parser.add_argument("--target-um", type=float, default=None,
                        help="Override target separation directly.")
    parser.add_argument("--export", choices=["env", "table"], default="table")
    args = parser.parse_args()

    # Determine target closest-point distance
    if args.target_um is None:
        bx, bz = SHAPES[args.baseline_shape]()
        target = closest_distance(bx, bz, args.pitch, args.baseline_gap)
    else:
        target = args.target_um

    if args.export == "table":
        print(f"# Pitch P = {args.pitch:g} µm,  target closest distance = {target:.3f} µm")
        print(f"# (matched to {args.baseline_shape} at initial_gap = {args.baseline_gap:g})")
        print(f"# DISP = |initial_gap| + 2*z_apex  (initial separation + full lobe diameter)")
        print(f"# {'shape':<12} {'initial_gap':>12} {'verified d':>12} {'z_apex':>8} {'DISP':>8}")
        for name, fn in SHAPES.items():
            hx, hz = fn()
            gap = find_gap_for_target(hx, hz, args.pitch, target)
            d = closest_distance(hx, hz, args.pitch, gap)
            z_apex = float(np.asarray(hz).max())
            disp = abs(gap) + 2.0 * z_apex
            print(f"  {name:<12} {gap:>12.3f} {d:>12.3f} {z_apex:>8.2f} {disp:>8.2f}")
    else:  # env
        for name, fn in SHAPES.items():
            hx, hz = fn()
            gap = find_gap_for_target(hx, hz, args.pitch, target)
            z_apex = float(np.asarray(hz).max())
            disp = abs(gap) + 2.0 * z_apex
            print(f"INITIAL_GAP_{name}={gap:.3f}")
            print(f"DISP_{name}={disp:.3f}")


if __name__ == "__main__":
    main()
