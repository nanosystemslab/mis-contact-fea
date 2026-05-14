# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "marimo",
#     "numpy",
#     "matplotlib",
#     "pyyaml",
# ]
# ///
"""Marimo + WASM port of the 2D alignment GUI.

Runs in two modes:

- **Local**: `marimo run src/mis_contact_fea/gui/align_marimo.py`
  (or `marimo edit` for the notebook UI).
- **Hosted**: GitHub Actions builds a static WASM bundle from this file
  and publishes it to GitHub Pages. The hosted page runs entirely in
  the user's browser via Pyodide — no server, no cold start, no idle
  sleep, no resource ceiling other than the user's own machine.

Self-contained design
---------------------
This notebook **inlines** the small subset of `mis_contact_fea` code it
needs (the seven shape generators and the polyline-layout helpers) so
it doesn't depend on the package being importable inside Pyodide. The
canonical source for these functions is:

    src/mis_contact_fea/profiles/builtins.py
    src/mis_contact_fea/mesh/layout.py

When you update those, **update this file too**. Tests in
`tests/test_marimo_inline_sync.py` enforce that the inlined SHAPES
dict produces byte-identical polylines to the package version.
"""
from __future__ import annotations

import marimo

__generated_with = "0.10.0"
app = marimo.App(width="medium", app_title="mis-contact-fea alignment")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        r"""
        # mis-contact-fea — 2D alignment

        Dial in your simulation visually, then **download the YAML** and run it locally:

        ```bash
        ./scripts/run_config.sh <downloaded>.yaml
        ```

        This page runs entirely in your browser (Pyodide WASM). Nothing is uploaded.
        """
    )
    return


@app.cell
def _():
    """Standard-library + numerics imports."""
    import io
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np
    import yaml
    return Path, io, np, plt, yaml


@app.cell
def _(np):
    """Inlined shape generators.

    SOURCE OF TRUTH: src/mis_contact_fea/profiles/builtins.py
    Keep in sync. The tests/test_marimo_inline_sync.py test pins
    byte-identical output to catch drift.
    """
    N_ARC = 80

    def gen_ellipse_axisym(Rxy, Rz, Rxy_p, Hp, n=N_ARC):
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

    def gen_cone_axisym(Rxy, H, Rxy_p, Hp, n_slope=N_ARC, fillet_r=0.0, n_fillet=20,
                        apex_fillet_r=0.0, n_apex_fillet=20):
        L = float(np.hypot(Rxy, H))
        slope_away = (Rxy / L, -H / L)
        beta = float(np.arctan2(Rxy, H))

        if apex_fillet_r > 0:
            r = apex_fillet_r
            apex_center = (0.0, H - r / np.sin(beta))
            t_along = r / np.tan(beta)
            apex_t_slope = (t_along * abs(slope_away[0]),
                            H + t_along * slope_away[1])
            apex_t_axis = (0.0, apex_center[1] + r)
            slope_end_x, slope_end_z = apex_t_slope
        else:
            slope_end_x, slope_end_z = 0.0, H

        if fillet_r <= 0:
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

        e2 = np.array([-Rxy / L, H / L])
        cos_alpha = -e2[0]
        alpha = float(np.arccos(cos_alpha))

        t_along = fillet_r / np.tan(alpha / 2.0)
        if t_along > (Rxy - Rxy_p) - 1e-9:
            raise ValueError(f"fillet_r={fillet_r} too large for cone")

        d = fillet_r / np.sin(alpha / 2.0)
        bis = np.array([-1.0, 0.0]) + e2
        bis = bis / float(np.linalg.norm(bis))
        cx, cz = (np.array([Rxy, 0.0]) + d * bis).tolist()

        t1 = np.array([cx, 0.0])
        t2 = np.array([Rxy + t_along * e2[0], t_along * e2[1]])

        pre_x = np.array([0.0, Rxy_p, Rxy_p, t1[0]])
        pre_z = np.array([-Hp,  -Hp,  0.0, 0.0])

        theta1 = float(np.arctan2(t1[1] - cz, t1[0] - cx))
        theta2 = float(np.arctan2(t2[1] - cz, t2[0] - cx))
        if theta2 < theta1:
            theta2 += 2.0 * np.pi
        thetas = np.linspace(theta1, theta2, n_fillet + 1)[1:]
        arc_x = cx + fillet_r * np.cos(thetas)
        arc_z = cz + fillet_r * np.sin(thetas)

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
        pillar_x = np.array([0.0, Rxy_p, Rxy_p])
        pillar_z = np.array([-Hp, -Hp, 0.0])
        t = np.linspace(0, np.pi / 2, n)
        arc_x = Rxy * np.cos(t)
        arc_z = Rz * np.sin(t)
        return np.concatenate([pillar_x, arc_x]), np.concatenate([pillar_z, arc_z])

    def gen_capcone_axisym(Rxy, Rz, r, Rxy_p, Hp, n=N_ARC):
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
        if Rxy_p >= Rxy:
            raise ValueError(f"Pillar must be narrower than doublecone: Rxy_p={Rxy_p}, Rxy={Rxy}")
        z_int = (H / 2) * (Rxy_p / Rxy - 1)

        if fillet_r <= 0:
            pillar_x = np.array([0.0, Rxy_p, Rxy_p])
            pillar_z = np.array([-Hp, -Hp, z_int])
            t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
            lower_x = Rxy_p + (Rxy - Rxy_p) * t
            lower_z = z_int + (0.0 - z_int) * t
            upper_x = Rxy * (1.0 - t)
            upper_z = (H / 2.0) * t
            return (np.concatenate([pillar_x, lower_x, upper_x]),
                    np.concatenate([pillar_z, lower_z, upper_z]))

        L_low = float(np.hypot(Rxy - Rxy_p, z_int))
        e_low = np.array([(Rxy_p - Rxy) / L_low, z_int / L_low])
        L_up = float(np.hypot(Rxy, H / 2.0))
        e_up = np.array([-Rxy / L_up, (H / 2.0) / L_up])

        cos_alpha = float(np.dot(e_low, e_up))
        alpha = float(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))

        t_along = fillet_r / np.tan(alpha / 2.0)
        if t_along > min(L_low, L_up) - 1e-9:
            raise ValueError(f"fillet_r={fillet_r} too large")

        d = fillet_r / np.sin(alpha / 2.0)
        bis = e_low + e_up
        bis = bis / float(np.linalg.norm(bis))
        cx, cz = (np.array([Rxy, 0.0]) + d * bis).tolist()

        t_lo = np.array([Rxy + t_along * e_low[0], t_along * e_low[1]])
        t_up = np.array([Rxy + t_along * e_up[0],  t_along * e_up[1]])

        pillar_x = [0.0, Rxy_p, Rxy_p]
        pillar_z = [-Hp, -Hp, z_int]

        t = np.linspace(0.0, 1.0, n_slope + 1)[1:]
        lower_x = Rxy_p + (t_lo[0] - Rxy_p) * t
        lower_z = z_int + (t_lo[1] - z_int) * t

        theta1 = float(np.arctan2(t_lo[1] - cz, t_lo[0] - cx))
        theta2 = float(np.arctan2(t_up[1] - cz, t_up[0] - cx))
        if theta2 < theta1:
            theta2 += 2.0 * np.pi
        thetas = np.linspace(theta1, theta2, n_fillet + 1)[1:]
        arc_x = cx + fillet_r * np.cos(thetas)
        arc_z = cz + fillet_r * np.sin(thetas)

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
        "sphere":     lambda Hp=220.0: gen_ellipse_axisym(50.0, 50.0, 25.0, Hp),
        "oblate":     lambda Hp=220.0: gen_ellipse_axisym(50.0, 25.0, 25.0, Hp),
        "prolate":    lambda Hp=220.0: gen_ellipse_axisym(50.0, 75.0, 25.0, Hp),
        "cone":       lambda Hp=220.0: gen_cone_axisym(53.033216, 50.0, 25.0, Hp, fillet_r=2.0, apex_fillet_r=2.0),
        "cap":        lambda Hp=220.0: gen_cap_axisym(50.0, 50.0, 25.0, Hp),
        "capcone":    lambda Hp=220.0: gen_capcone_axisym(50.0, 50.0, 25.0, 25.0, Hp),
        "doublecone": lambda Hp=220.0: gen_doublecone_axisym(50.848528, 100.0, 25.0, Hp, fillet_r=2.0, apex_fillet_r=2.0),
    }
    return (SHAPES,)


@app.cell
def _(np):
    """Inlined polyline-layout helpers.

    SOURCE OF TRUTH: src/mis_contact_fea/mesh/layout.py.
    """
    def bottom_body_vertices(half_x, half_z, P, plate_thickness):
        Hp = -float(half_z[0])
        plate_bot = -Hp - plate_thickness
        plate_top = -Hp
        n_poly = len(half_x)
        vx, vz, tags = [], [], []
        def add(x, z, edge_tag):
            vx.append(float(x)); vz.append(float(z)); tags.append(edge_tag)
        add(0.0,        plate_bot, "bottom_fixed")
        add(P,          plate_bot, "sym_right")
        add(P,          plate_top, "sym_right")
        add(P,          float(half_z[-1]), "bottom_contact")
        mir_x = (P - np.asarray(half_x))[::-1]
        mir_z = np.asarray(half_z)[::-1]
        for i in range(1, n_poly - 1):
            add(mir_x[i], mir_z[i], "bottom_contact")
        for i in range(1, n_poly - 1):
            add(float(half_x[i]), float(half_z[i]), "bottom_contact")
        add(float(half_x[-1]), float(half_z[-1]), "sym_left")
        add(0.0, plate_top, "sym_left")
        return np.array(vx), np.array(vz), tags

    def top_body_vertices(half_x, half_z, P, plate_thickness, initial_gap):
        Hp = -float(half_z[0])
        z_apex = float(np.asarray(half_z).max())
        inner = Hp + 2.0 * z_apex + initial_gap
        outer = inner + plate_thickness
        top_x_right = P / 2.0 + np.asarray(half_x)
        top_z = 2.0 * z_apex + initial_gap - np.asarray(half_z)
        n_poly = len(half_x)
        vx, vz, tags = [], [], []
        def add(x, z, edge_tag):
            vx.append(float(x)); vz.append(float(z)); tags.append(edge_tag)
        add(0.0, outer, "top_disp")
        add(P,   outer, "sym_right")
        add(P,   inner, "top_contact")
        for i in range(1, n_poly):
            add(float(top_x_right[i]), float(top_z[i]), "top_contact")
        top_x_left = P / 2.0 - np.asarray(half_x)
        for i in range(n_poly - 2, 0, -1):
            add(float(top_x_left[i]), float(top_z[i]), "top_contact")
        add(0.0, inner, "sym_left")
        return np.array(vx), np.array(vz), tags

    return bottom_body_vertices, top_body_vertices


@app.cell
def _():
    """Tuned default INITIAL_GAP / DISP / STEPS per shape — matches
    examples/01_seven_builtins."""
    RECOMMENDED = {
        #            (initial_gap, push_disp, push_steps)
        "sphere":     (-55.0,  85.0, 425),
        "oblate":     (-30.0,  40.0, 200),
        "prolate":    (-80.0,  87.0, 435),
        "cone":       (-75.0,  30.9, 155),
        "cap":        (-65.0,  40.0, 200),
        "capcone":   (-110.0,  45.0, 225),
        "doublecone": (-80.0,  40.0, 200),
    }
    return (RECOMMENDED,)


@app.cell
def _(SHAPES, mo):
    shape = mo.ui.dropdown(options=sorted(SHAPES.keys()), value="sphere", label="Shape")
    direction = mo.ui.radio(
        options=["down (push-in)", "up (retention)"],
        value="down (push-in)",
        label="Direction",
    )
    return direction, shape


@app.cell
def _(RECOMMENDED, mo, shape):
    rec_ig, rec_disp, rec_steps = RECOMMENDED[shape.value]
    initial_gap = mo.ui.slider(
        start=-250.0, stop=50.0, step=0.5, value=float(rec_ig),
        label="Initial gap (µm)  —  negative = bodies overlap at start",
    )
    disp = mo.ui.slider(
        start=0.5, stop=300.0, step=0.5, value=float(rec_disp),
        label="Total displacement (µm)",
    )
    steps = mo.ui.number(
        start=1, stop=5000, step=5, value=int(rec_steps),
        label="Newton increments (~0.2 µm/step recommended)",
    )
    return disp, initial_gap, steps


@app.cell
def _(direction, disp, initial_gap, mo, shape, steps):
    mo.vstack([
        mo.hstack([shape, direction]),
        initial_gap,
        disp,
        steps,
    ])
    return


@app.cell
def _(mo):
    pitch = mo.ui.slider(start=50.0, stop=400.0, step=1.0, value=198.0, label="Pitch (µm)")
    Hp = mo.ui.slider(start=80.0, stop=400.0, step=5.0, value=220.0, label="Pillar height Hp (µm)")
    plate = mo.ui.slider(start=5.0, stop=100.0, step=1.0, value=30.0, label="Plate thickness (µm)")
    mesh_h = mo.ui.slider(start=1.0, stop=20.0, step=0.5, value=4.0, label="Mesh size (µm)")
    return Hp, mesh_h, pitch, plate


@app.cell
def _(Hp, mesh_h, mo, pitch, plate):
    mo.accordion({
        "Geometry (advanced)": mo.vstack([pitch, Hp, plate, mesh_h]),
    })
    return


@app.cell
def _(mo):
    contact_mode = mo.ui.dropdown(
        options=["raytracing", "closest"], value="raytracing",
        label="Contact mode (raytracing for push, closest for retention)",
    )
    fric = mo.ui.number(start=0.0, stop=2.0, step=0.05, value=0.0, label="Friction coefficient")
    newton_rtol = mo.ui.number(start=1e-6, stop=1.0, step=0.01, value=0.2, label="Newton rtol")
    newton_atol = mo.ui.number(start=1e-6, stop=1e6, step=10.0, value=1000.0, label="Newton atol")
    gamma_scale = mo.ui.number(start=0.001, stop=100.0, step=0.05, value=0.2, label="Gamma scale (Nitsche penalty)")
    return contact_mode, fric, gamma_scale, newton_atol, newton_rtol


@app.cell
def _(contact_mode, fric, gamma_scale, mo, newton_atol, newton_rtol):
    mo.accordion({
        "Solver (advanced)": mo.vstack([contact_mode, fric, newton_rtol, newton_atol, gamma_scale]),
    })
    return


@app.cell
def _(direction, mo, shape):
    dir_value = "down" if direction.value.startswith("down") else "up"
    default_out_dir = (
        f"runs_all_shapes/{shape.value}" if dir_value == "down"
        else f"runs_retention/{shape.value}"
    )
    out_dir = mo.ui.text(value=default_out_dir, label="Output directory")
    return dir_value, out_dir


@app.cell
def _(out_dir):
    out_dir
    return


@app.cell
def _(
    Hp,
    SHAPES,
    bottom_body_vertices,
    dir_value,
    disp,
    initial_gap,
    np,
    pitch,
    plate,
    plt,
    shape,
    top_body_vertices,
):
    """Two-panel matplotlib preview."""
    half_x, half_z = SHAPES[shape.value](Hp=Hp.value)
    bvx, bvz, _ = bottom_body_vertices(half_x, half_z, pitch.value, plate.value)
    tvx, tvz_init, _ = top_body_vertices(half_x, half_z, pitch.value, plate.value, initial_gap.value)
    bvx, bvz = np.asarray(bvx), np.asarray(bvz)
    tvx, tvz_init = np.asarray(tvx), np.asarray(tvz_init)

    sign = -1.0 if dir_value == "down" else +1.0
    tvz_final = tvz_init + sign * disp.value

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=110)
    z_lo = min(bvz.min(), tvz_init.min(), tvz_final.min()) - 10
    z_hi = max(bvz.max(), tvz_init.max(), tvz_final.max()) + 10

    for ax, tvz, title in [
        (axes[0], tvz_init,  "initial   (top plate disp = 0)"),
        (axes[1], tvz_final, f"final     (top plate disp = {sign * disp.value:+.2f} µm)"),
    ]:
        ax.fill(bvx, bvz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.2)
        ax.fill(tvx, tvz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.2)
        ax.axvline(0,            linestyle="--", color="#888", linewidth=0.8)
        ax.axvline(pitch.value,  linestyle="--", color="#888", linewidth=0.8)
        ax.set_aspect("equal")
        ax.set_xlim(-15, pitch.value + 15)
        ax.set_ylim(z_lo, z_hi)
        ax.grid(alpha=0.2)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("z (µm)")
        ax.set_title(title, fontsize=10)
    fig.suptitle(
        f"{shape.value} — IG = {initial_gap.value:+g} µm, DISP = {disp.value:g} µm, pitch = {pitch.value:g} µm",
        fontsize=10,
    )
    fig.tight_layout()
    fig
    return


@app.cell
def _(
    Hp,
    contact_mode,
    dir_value,
    disp,
    fric,
    gamma_scale,
    initial_gap,
    io,
    mesh_h,
    newton_atol,
    newton_rtol,
    out_dir,
    pitch,
    plate,
    shape,
    steps,
    yaml,
):
    """Build the YAML config from current UI state (plain dict — no pydantic
    needed; the local CLI revalidates with the full pydantic schema)."""
    cfg_dict = {
        "geometry": {
            "pitch_um": float(pitch.value),
            "Hp_um": float(Hp.value),
            "plate_um": float(plate.value),
            "initial_gap_um": float(initial_gap.value),
        },
        "profile": {"kind": "builtin", "name": shape.value, "units": "um"},
        "mesh": {"characteristic_size_um": float(mesh_h.value)},
        "solver": {
            "steps": int(steps.value),
            "disp_um": float(disp.value),
            "direction": dir_value,
            "contact_mode": contact_mode.value,
            "fric": float(fric.value),
            "newton_rtol": float(newton_rtol.value),
            "newton_atol": float(newton_atol.value),
            "gamma_scale": float(gamma_scale.value),
            "ksp_type": "preonly",
            "pc_type": "lu",
        },
        "output": {"out_dir": out_dir.value},
    }
    buf = io.StringIO()
    yaml.safe_dump(cfg_dict, buf, sort_keys=False, default_flow_style=False)
    yaml_text = buf.getvalue()
    return (yaml_text,)


@app.cell
def _(mo, yaml_text):
    mo.md(f"### YAML preview\n\n```yaml\n{yaml_text}\n```")
    return


@app.cell
def _(dir_value, mo, shape, yaml_text):
    fname = f"{shape.value}_{dir_value}.yaml"
    download = mo.download(
        data=yaml_text.encode("utf-8"),
        filename=fname,
        mimetype="text/yaml",
        label=f"Download {fname}",
    )
    download
    return (fname,)


@app.cell
def _(fname, mo):
    mo.md(
        rf"""
        ---
        ### Running the downloaded config

        Save `{fname}` somewhere in your local checkout of the
        `mis-contact-fea` repo, then:

        ```bash
        ./scripts/run_config.sh path/to/{fname}
        ```

        Or batch several at once:

        ```bash
        ./scripts/run_configs.sh path/to/*.yaml
        ```

        See [examples/](https://github.com/nanosystemslab/mis-contact-fea/tree/main/examples)
        for the full sweep recipes.
        """
    )
    return


if __name__ == "__main__":
    app.run()
