"""Right-half axisymmetric profile generators, lifted from
family_2d_geometries_axisym.ipynb so they can be imported by the 2D
plane-strain mesher.

Each generator returns (x, z) of a polyline starting at (0, z_bottom) on
the axis and ending at (0, z_top) on the axis. The polyline traces the
right exterior of (pillar ∪ structure):
    index 0: (0, -Hp)        axis-base of pillar
    index 1: (Rxy_p, -Hp)    pillar bottom right corner
    index 2: (Rxy_p, z_int)  pillar shaft transitions into structure
    ...
    index N: (0, z_top)      structure apex on axis

The bottom-edge segment (index 0 -> index 1) is interior to the merged
body once a plate is glued below, so the slice mesher skips index 0
when building outlines.
"""
from __future__ import annotations

import numpy as np

N_ARC = 80


def gen_ellipse_axisym(Rxy, Rz, Rxy_p, Hp, n=N_ARC):
    """Right-half ellipsoid centered at (0, 0) merged with pillar below.
    Pillar: x in [0, Rxy_p], z in [-Hp, 0]. The pillar's upper portion
    is volumetrically inside the ellipsoid's lower half and is hidden
    from the merged outline."""
    if Rxy_p >= Rxy:
        raise ValueError(f"Pillar must be narrower than ellipse: Rxy_p={Rxy_p}, Rxy={Rxy}")
    z_int = -Rz * np.sqrt(1 - (Rxy_p / Rxy) ** 2)
    pillar_x = np.array([0.0, Rxy_p])
    pillar_z = np.array([-Hp, -Hp])
    theta_start = np.arctan2(z_int / Rz, Rxy_p / Rxy)
    theta_end = np.pi / 2
    t = np.linspace(theta_start, theta_end, n)
    arc_x = Rxy * np.cos(t)
    arc_z = Rz * np.sin(t)
    return np.concatenate([pillar_x, arc_x]), np.concatenate([pillar_z, arc_z])


def _apex_fillet_arc(slope_end_dir, apex_xy, r, n):
    """Return arc points (x, z) for a fillet at an axis-apex (right-half).

    The apex is at apex_xy = (0, z_apex). slope_end_dir is the unit vector
    along the slope direction AWAY from the apex (i.e., toward the base).
    The fillet replaces the sharp apex with a quarter-ish arc tangent to
    the slope on one end and to the axis (x=0 line) on the other.

    The full cross-section's interior angle at the apex is 2*beta where
    beta = angle between slope and axis = arctan(slope_end_dir.x /
    |slope_end_dir.z|). For our symmetric cone/doublecone shapes the
    apex angle is 2*beta and the fillet center sits on the axis at
    z = z_apex - r/sin(beta).
    """
    apex_z = float(apex_xy[1])
    # slope_end_dir points from apex toward base; its z-component is
    # negative (going downward). beta is the angle between slope and axis.
    dx, dz = float(slope_end_dir[0]), float(slope_end_dir[1])
    beta = np.arctan2(abs(dx), abs(dz))
    if beta <= 1e-9:
        return np.array([0.0]), np.array([apex_z - r])
    # Tangent length along slope from apex (= r / tan(beta) for symmetric apex)
    # Center of fillet on the axis, below apex by r / sin(beta).
    center = np.array([0.0, apex_z - r / np.sin(beta)])
    # Slope tangent point in the right half:
    t_slope = np.array([0.0 + (r / np.tan(beta)) * (dx / abs(dx) if dx != 0 else 0.0),
                        apex_z + (r / np.tan(beta)) * dz / abs(dz) * (-1.0 if dz < 0 else 1.0)])
    # Cleaner: tangent point = apex + (r/tan(beta)) * slope_end_dir
    t_slope = np.array([apex_xy[0] + (r / np.tan(beta)) * dx,
                        apex_xy[1] + (r / np.tan(beta)) * dz])
    # Axis tangent point: on the axis (x=0) at z = center.z + r
    t_axis = np.array([0.0, center[1] + r])
    # Arc from t_slope to t_axis (CCW, with center as origin).
    theta1 = float(np.arctan2(t_slope[1] - center[1], t_slope[0] - center[0]))
    theta2 = float(np.arctan2(t_axis[1] - center[1], t_axis[0] - center[0]))
    if theta2 < theta1:
        theta2 += 2.0 * np.pi
    thetas = np.linspace(theta1, theta2, n + 1)[1:]  # skip start (= t_slope)
    return (center[0] + r * np.cos(thetas), center[1] + r * np.sin(thetas))


def gen_cone_axisym(Rxy, H, Rxy_p, Hp, n_slope=N_ARC, fillet_r=0.0, n_fillet=20,
                    apex_fillet_r=0.0, n_apex_fillet=20):
    """Right-half cone (base at z=0, apex at (0, H)) + pillar below.

    Optionally insert a circular fillet at the base/slope corner
    (Rxy, 0) and at the apex (0, H). Both points are stress
    concentrations for sharp linear-elastic contact; rounding them
    makes Newton convergence stable.

    Parameters
    ----------
    Rxy, H : outer base radius and apex height of the un-filleted cone
    Rxy_p, Hp : pillar shaft radius and height (z in [-Hp, 0])
    n_slope : densification of the cone slope segment
    fillet_r : fillet radius at (Rxy, 0) base/slope corner. 0 = sharp.
    n_fillet : number of points along the base fillet arc
    apex_fillet_r : fillet radius at the apex (0, H). 0 = sharp tip.
    n_apex_fillet : number of points along the apex fillet arc
    """
    # Compute slope-direction unit vector AWAY from apex (toward base).
    L = float(np.hypot(Rxy, H))
    slope_away = (Rxy / L, -H / L)
    # beta = angle between slope and axis (= half of the apex's interior
    # angle in the full cross-section). For symmetric cone:
    beta = float(np.arctan2(Rxy, H))  # arctan(Rxy/H)

    # Apex fillet: a small hemisphere centered on the axis, tangent to
    # the cone slope. This is the geometrically natural "rounded cone
    # tip". There IS a 90° corner where the arc meets the cell-edge at
    # the axis (arc tangent there is horizontal, cell-edge is vertical),
    # but trying to make the arc smooth-tangent to the axis creates an
    # unphysical bulge that's worse for stress.
    if apex_fillet_r > 0:
        r = apex_fillet_r
        # Center on axis, at distance r/sin(beta) below the apex.
        apex_center = (0.0, H - r / np.sin(beta))
        # Slope tangent: along the slope, at distance r/tan(beta) from apex.
        t_along = r / np.tan(beta)
        apex_t_slope = (t_along * abs(slope_away[0]),
                        H + t_along * slope_away[1])
        apex_t_axis = (0.0, apex_center[1] + r)
        slope_end_x, slope_end_z = apex_t_slope
    else:
        slope_end_x, slope_end_z = 0.0, H

    if fillet_r <= 0:
        # Sharp base corner.
        pillar_x = np.array([0.0, Rxy_p, Rxy_p, Rxy])
        pillar_z = np.array([-Hp, -Hp, 0.0, 0.0])
        slope_start_x, slope_start_z = Rxy, 0.0
        t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
        slope_x = slope_start_x + (slope_end_x - slope_start_x) * t
        slope_z = slope_start_z + (slope_end_z - slope_start_z) * t
        out_x = [pillar_x, slope_x]
        out_z = [pillar_z, slope_z]
        if apex_fillet_r > 0:
            thetas = np.linspace(
                np.arctan2(apex_t_slope[1] - apex_center[1], apex_t_slope[0] - apex_center[0]),
                np.arctan2(apex_t_axis[1]  - apex_center[1], apex_t_axis[0]  - apex_center[0]),
                n_apex_fillet + 1,
            )
            if thetas[-1] < thetas[0]:
                thetas[1:] += 2.0 * np.pi
            out_x.append(apex_center[0] + apex_fillet_r * np.cos(thetas[1:]))
            out_z.append(apex_center[1] + apex_fillet_r * np.sin(thetas[1:]))
        return np.concatenate(out_x), np.concatenate(out_z)

    # Geometry at the base/slope corner (Rxy, 0):
    #   Edge 1 (horizontal toward pillar): direction (-1, 0)
    #   Edge 2 (slope toward apex):       direction (-Rxy, H)/L
    L = float(np.hypot(Rxy, H))
    e2 = np.array([-Rxy / L, H / L])
    cos_alpha = -e2[0]  # = Rxy/L, between 0 and 1
    alpha = float(np.arccos(cos_alpha))  # interior angle at the corner

    # Sanity: fillet must fit between (Rxy_p, 0) and the apex along its tangents.
    t_along = fillet_r / np.tan(alpha / 2.0)
    if t_along > (Rxy - Rxy_p) - 1e-9:
        raise ValueError(
            f"fillet_r={fillet_r} too large for cone (Rxy={Rxy}, Rxy_p={Rxy_p}, H={H}); "
            f"tangent extent {t_along:.2f} exceeds horizontal length {Rxy - Rxy_p:.2f}."
        )

    # Center of the fillet on the angle bisector, distance d from the corner.
    d = fillet_r / np.sin(alpha / 2.0)
    bis = np.array([-1.0, 0.0]) + e2
    bis = bis / float(np.linalg.norm(bis))
    cx, cz = (np.array([Rxy, 0.0]) + d * bis).tolist()

    # Tangent points on each edge.
    t1 = np.array([cx, 0.0])                                      # on horizontal
    t2 = np.array([Rxy + t_along * e2[0], t_along * e2[1]])       # on slope

    # Polyline assembly:
    #   axis-base, pillar-base-right, pillar-shaft-top, horizontal_end (=t1),
    #   fillet arc points (skip first = t1), slope densification ending at apex.
    pre_x = np.array([0.0, Rxy_p, Rxy_p, t1[0]])
    pre_z = np.array([-Hp,  -Hp,  0.0, 0.0])

    # Fillet arc from t1 to t2 (CCW since the arc bulges INTO the body).
    theta1 = float(np.arctan2(t1[1] - cz, t1[0] - cx))
    theta2 = float(np.arctan2(t2[1] - cz, t2[0] - cx))
    if theta2 < theta1:
        theta2 += 2.0 * np.pi
    thetas = np.linspace(theta1, theta2, n_fillet + 1)[1:]  # skip start (= t1)
    arc_x = cx + fillet_r * np.cos(thetas)
    arc_z = cz + fillet_r * np.sin(thetas)

    # Slope densified from t2 to slope_end (apex or apex-fillet tangent).
    t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
    slope_x = t2[0] + (slope_end_x - t2[0]) * t
    slope_z = t2[1] + (slope_end_z - t2[1]) * t

    out_x = [pre_x, arc_x, slope_x]
    out_z = [pre_z, arc_z, slope_z]

    if apex_fillet_r > 0:
        thetas = np.linspace(
            np.arctan2(apex_t_slope[1] - apex_center[1], apex_t_slope[0] - apex_center[0]),
            np.arctan2(apex_t_axis[1]  - apex_center[1], apex_t_axis[0]  - apex_center[0]),
            n_apex_fillet + 1,
        )
        if thetas[-1] < thetas[0]:
            thetas[1:] += 2.0 * np.pi
        out_x.append(apex_center[0] + apex_fillet_r * np.cos(thetas[1:]))
        out_z.append(apex_center[1] + apex_fillet_r * np.sin(thetas[1:]))

    return np.concatenate(out_x), np.concatenate(out_z)


def gen_cap_axisym(Rxy, Rz, Rxy_p, Hp, n=N_ARC):
    """Right-half hemi-ellipse cap (base at z=0, apex at (0, Rz)) + pillar.
    Flat-bottomed: cap and pillar meet at z=0."""
    pillar_x = np.array([0.0, Rxy_p, Rxy_p])
    pillar_z = np.array([-Hp, -Hp, 0.0])
    t = np.linspace(0, np.pi / 2, n)
    arc_x = Rxy * np.cos(t)
    arc_z = Rz * np.sin(t)
    return np.concatenate([pillar_x, arc_x]), np.concatenate([pillar_z, arc_z])


def gen_capcone_axisym(Rxy, Rz, r, Rxy_p, Hp, n=N_ARC):
    """Right-half cap-cone (umbonate hemisphere-compound) + pillar below.
    Currently supports only Rxy == Rz."""
    if Rxy != Rz:
        raise NotImplementedError("gen_capcone_axisym currently supports only Rxy == Rz")
    R_main, R_small, d = Rxy, r, Rz
    sin_theta = (R_main - R_small) / d
    theta = np.arcsin(sin_theta)
    pillar_x = np.array([0.0, Rxy_p, Rxy_p])
    pillar_z = np.array([-Hp, -Hp, 0.0])
    n_main = max(int(n * theta / (np.pi / 2)), 5)
    t_main = np.linspace(0, theta, n_main)
    main_x = R_main * np.cos(t_main)
    main_z = R_main * np.sin(t_main)
    n_small = max(int(n * (np.pi / 2 - theta) / (np.pi / 2)), 10)
    t_small = np.linspace(theta, np.pi / 2, n_small)
    small_x = R_small * np.cos(t_small)
    small_z = d + R_small * np.sin(t_small)
    x = np.concatenate([pillar_x, main_x, small_x])
    z = np.concatenate([pillar_z, main_z, small_z])
    return x, z


def gen_doublecone_axisym(Rxy, H, Rxy_p, Hp, n_slope=N_ARC, fillet_r=0.0, n_fillet=20,
                          apex_fillet_r=0.0, n_apex_fillet=20):
    """Right-half doublecone (rhombus) centered at (0, 0) + pillar below.

    Optionally insert a circular fillet at the equator corner (Rxy, 0)
    where the lower and upper slopes meet. Interior angle there is 90°,
    and rounding it eliminates the dominant stress-concentration source
    just as the cone's base/slope fillet does.

    Both slopes are densified so closest-point projection has well-
    defined unique pairs."""
    if Rxy_p >= Rxy:
        raise ValueError(f"Pillar must be narrower than doublecone: Rxy_p={Rxy_p}, Rxy={Rxy}")
    z_int = (H / 2) * (Rxy_p / Rxy - 1)

    if fillet_r <= 0:
        # Sharp equator corner — legacy behavior.
        pillar_x = np.array([0.0, Rxy_p, Rxy_p])
        pillar_z = np.array([-Hp, -Hp, z_int])
        t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
        lower_x = Rxy_p + (Rxy - Rxy_p) * t
        lower_z = z_int + (0.0 - z_int) * t
        upper_x = Rxy * (1.0 - t)
        upper_z = (H / 2.0) * t
        return (
            np.concatenate([pillar_x, lower_x, upper_x]),
            np.concatenate([pillar_z, lower_z, upper_z]),
        )

    # Geometry at equator corner (Rxy, 0):
    #   Lower slope direction away from corner: toward (Rxy_p, z_int)
    #   Upper slope direction away from corner: toward (0, H/2)
    L_low = float(np.hypot(Rxy - Rxy_p, z_int))
    e_low = np.array([(Rxy_p - Rxy) / L_low, z_int / L_low])
    L_up = float(np.hypot(Rxy, H / 2.0))
    e_up = np.array([-Rxy / L_up, (H / 2.0) / L_up])

    cos_alpha = float(np.dot(e_low, e_up))
    alpha = float(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))  # interior angle

    t_along = fillet_r / np.tan(alpha / 2.0)
    if t_along > min(L_low, L_up) - 1e-9:
        raise ValueError(
            f"fillet_r={fillet_r} too large for doublecone geometry "
            f"(slope lengths {L_low:.2f}, {L_up:.2f}); tangent extent {t_along:.2f}."
        )

    d = fillet_r / np.sin(alpha / 2.0)
    bis = e_low + e_up
    bis = bis / float(np.linalg.norm(bis))
    cx, cz = (np.array([Rxy, 0.0]) + d * bis).tolist()

    t_lo = np.array([Rxy + t_along * e_low[0], t_along * e_low[1]])  # on lower slope
    t_up = np.array([Rxy + t_along * e_up[0],  t_along * e_up[1]])   # on upper slope

    # Pillar shaft up to lower slope start (unchanged geometry below the fillet).
    pillar_x = [0.0, Rxy_p, Rxy_p]
    pillar_z = [-Hp, -Hp, z_int]

    # Lower slope densified from (Rxy_p, z_int) to t_lo (just below equator).
    t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
    lower_x = Rxy_p + (t_lo[0] - Rxy_p) * t
    lower_z = z_int + (t_lo[1] - z_int) * t

    # Fillet arc t_lo -> t_up, CCW (bulges into body which is to the left).
    theta1 = float(np.arctan2(t_lo[1] - cz, t_lo[0] - cx))
    theta2 = float(np.arctan2(t_up[1] - cz, t_up[0] - cx))
    if theta2 < theta1:
        theta2 += 2.0 * np.pi
    thetas = np.linspace(theta1, theta2, n_fillet + 1)[1:]
    arc_x = cx + fillet_r * np.cos(thetas)
    arc_z = cz + fillet_r * np.sin(thetas)

    # Apex fillet at (0, H/2): hemisphere centered on axis, tangent to
    # the upper slope. There IS a 90° corner where the arc meets the
    # cell-edge at the axis, but this is the natural rounded-tip shape.
    H_half = H / 2.0
    Lup = float(np.hypot(Rxy, H_half))
    beta_apex_dc = float(np.arctan2(Rxy, H_half))
    if apex_fillet_r > 0:
        r = apex_fillet_r
        apex_center = (0.0, H_half - r / np.sin(beta_apex_dc))
        t_along = r / np.tan(beta_apex_dc)
        apex_t_slope = (t_along * (Rxy / Lup),
                        H_half + t_along * (-H_half / Lup))
        apex_t_axis = (0.0, apex_center[1] + r)
        upper_end_x, upper_end_z = apex_t_slope
    else:
        upper_end_x, upper_end_z = 0.0, H_half

    # Upper slope densified from t_up to upper_end (apex or apex-fillet tangent).
    upper_x = t_up[0] + (upper_end_x - t_up[0]) * t
    upper_z = t_up[1] + (upper_end_z - t_up[1]) * t

    out_x = [pillar_x, lower_x, arc_x, upper_x]
    out_z = [pillar_z, lower_z, arc_z, upper_z]

    if apex_fillet_r > 0:
        thetas = np.linspace(
            np.arctan2(apex_t_slope[1] - apex_center[1], apex_t_slope[0] - apex_center[0]),
            np.arctan2(apex_t_axis[1]  - apex_center[1], apex_t_axis[0]  - apex_center[0]),
            n_apex_fillet + 1,
        )
        if thetas[-1] < thetas[0]:
            thetas[1:] += 2.0 * np.pi
        out_x.append(apex_center[0] + apex_fillet_r * np.cos(thetas[1:]))
        out_z.append(apex_center[1] + apex_fillet_r * np.sin(thetas[1:]))

    return np.concatenate(out_x), np.concatenate(out_z)


SHAPES = {
    "sphere":      lambda Hp=220.0: gen_ellipse_axisym(50.0, 50.0, 25.0, Hp),
    "oblate":      lambda Hp=220.0: gen_ellipse_axisym(50.0, 25.0, 25.0, Hp),
    "prolate":     lambda Hp=220.0: gen_ellipse_axisym(50.0, 75.0, 25.0, Hp),
    # Cone/doublecone construction Rxy is shrunk so the OUTER widest point of the
    # filleted body matches Rxy=50 for the ellipsoid family. With fillet_r=2 µm,
    # the rounded corner bulges out of the construction radius; the offsets below
    # were computed by bisection to give max_x = 49.996 µm (same as ellipses).
    "cone":        lambda Hp=220.0: gen_cone_axisym(53.033216, 50.0, 25.0, Hp, fillet_r=2.0, apex_fillet_r=2.0),
    "cap":         lambda Hp=220.0: gen_cap_axisym(50.0, 50.0, 25.0, Hp),
    "capcone":     lambda Hp=220.0: gen_capcone_axisym(50.0, 50.0, 25.0, 25.0, Hp),
    "doublecone":  lambda Hp=220.0: gen_doublecone_axisym(50.848528, 100.0, 25.0, Hp, fillet_r=2.0, apex_fillet_r=2.0),
}
