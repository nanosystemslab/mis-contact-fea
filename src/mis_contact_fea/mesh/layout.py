"""Pure-numpy polyline-layout helpers for the 2D unit cell.

These functions compute the closed-polygon outlines of the bottom and
top bodies given a right-half axisymmetric profile, plus a few
geometric utilities. They have **no gmsh dependency** so they can be
imported in environments where gmsh isn't available — notably the
WASM-hosted Marimo GUI and the Streamlit alignment app.

The actual mesh assembly (gmsh) lives next door in `slice_2d.py`.
"""
from __future__ import annotations

import numpy as np


def bottom_body_vertices(half_x, half_z, P, plate_thickness):
    """CCW vertices of the bottom body's outline as parallel arrays.

    Also returns a list `edge_tags` of length len(verts), where
    edge_tags[i] is the physical-group name for the edge from vert i
    to vert (i+1) % len(verts).
    """
    Hp = -float(half_z[0])
    plate_bot = -Hp - plate_thickness
    plate_top = -Hp

    # half_x has N+1 entries (indices 0..N). Original polyline starts at
    # (0, plate_top), ends at (0, z_apex).
    # Mirrored polyline = (P - half_x, half_z), same N+1 entries.
    n_poly = len(half_x)  # = N + 1

    vx, vz, tags = [], [], []

    def add(x, z, edge_tag):
        vx.append(float(x))
        vz.append(float(z))
        tags.append(edge_tag)

    # Plate corners + cell-right-edge cut
    add(0.0,        plate_bot, "bottom_fixed")  # V0
    add(P,          plate_bot, "sym_right")     # V1
    add(P,          plate_top, "sym_right")     # V2
    add(P,          float(half_z[-1]), "bottom_contact")  # V3 = (P, z_apex)

    # Mirrored polyline (curve of pillar half-B's left side, traversed
    # top->bottom).
    mir_x = (P - np.asarray(half_x))[::-1]
    mir_z = np.asarray(half_z)[::-1]
    for i in range(1, n_poly - 1):
        add(mir_x[i], mir_z[i], "bottom_contact")

    # Now traverse the original polyline indices 1..N-1 (pillar half-A
    # right side, bottom->top, excluding the axis-base and the apex).
    for i in range(1, n_poly - 1):
        add(float(half_x[i]), float(half_z[i]), "bottom_contact")
    # Last polyline point = apex on axis, also where sym_left starts.
    add(float(half_x[-1]), float(half_z[-1]), "sym_left")

    # Intermediate vertex at (0, plate_top) to break the long sym_left
    # edge into two segments.
    add(0.0, plate_top, "sym_left")

    return np.array(vx), np.array(vz), tags


def top_body_vertices(half_x, half_z, P, plate_thickness, initial_gap):
    """CCW vertices of the top body's outline + edge tags."""
    Hp = -float(half_z[0])
    z_apex = float(np.asarray(half_z).max())
    inner = Hp + 2.0 * z_apex + initial_gap  # plate-bottom z (plate sits above)
    outer = inner + plate_thickness

    # Mirrored-in-z polyline of top pillar (in cell coords, centered at x=P/2)
    top_x_right = P / 2.0 + np.asarray(half_x)  # right side of top pillar
    top_z = 2.0 * z_apex + initial_gap - np.asarray(half_z)
    n_poly = len(half_x)

    vx, vz, tags = [], [], []

    def add(x, z, edge_tag):
        vx.append(float(x))
        vz.append(float(z))
        tags.append(edge_tag)

    add(0.0, outer, "top_disp")    # V0
    add(P,   outer, "sym_right")   # V1
    add(P,   inner, "top_contact")
    for i in range(1, n_poly):
        add(float(top_x_right[i]), float(top_z[i]), "top_contact")
    top_x_left = P / 2.0 - np.asarray(half_x)
    for i in range(n_poly - 2, 0, -1):
        add(float(top_x_left[i]), float(top_z[i]), "top_contact")
    add(0.0, inner, "sym_left")

    return np.array(vx), np.array(vz), tags


def signed_area(vx, vz) -> float:
    """Signed area of a closed polygon; positive = CCW."""
    n = len(vx)
    s = 0.0
    for i in range(n):
        j = (i + 1) % n
        s += vx[i] * vz[j] - vx[j] * vz[i]
    return 0.5 * s


def dedup_consecutive(vx, vz, edge_tags, tol: float = 1e-9):
    """Drop any vertex coincident with its successor (within tol).

    A duplicate creates a zero-length edge in the Gmsh geometry, which
    causes the 2D mesher to hang. The pillar/arc-transition in some
    profiles introduces such duplicates.
    """
    n = len(vx)
    keep = []
    for i in range(n):
        j = (i + 1) % n
        if abs(vx[i] - vx[j]) < tol and abs(vz[i] - vz[j]) < tol:
            continue  # drop vertex i; its outgoing edge is zero-length
        keep.append(i)
    return (
        [vx[i] for i in keep],
        [vz[i] for i in keep],
        [edge_tags[i] for i in keep],
    )
