# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "marimo",
#     "numpy",
#     "matplotlib",
#     "pyyaml",
#     "svgpathtools",
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
    """Inlined bytes-based loaders for CSV / SVG profile uploads.

    The package versions in `mis_contact_fea.profiles.from_csv` and
    `from_svg` take filesystem paths; in the browser we receive bytes
    from `mo.ui.file`, so we re-implement against `io.StringIO`. Mirrors
    the same polyline contract (start + end on the axis, x >= 0, µm).
    """
    import csv as _csv
    import io as _io

    _UM_PER_UNIT = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1_000_000.0}

    def _looks_like_header(row):
        try:
            float(row[0])
            return False
        except (ValueError, IndexError):
            return True

    def load_csv_bytes(data, Hp=220.0, units="um", prepend_pillar=False, Rxy_p=25.0):
        scale = _UM_PER_UNIT[units]
        text = data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else data
        reader = _csv.reader(_io.StringIO(text))
        rows_raw = [r for r in reader if r and not r[0].startswith("#")]
        if not rows_raw:
            raise ValueError("CSV is empty")
        if _looks_like_header(rows_raw[0]):
            rows_raw = rows_raw[1:]
        x = np.array([float(r[0]) for r in rows_raw]) * scale
        z = np.array([float(r[1]) for r in rows_raw]) * scale
        if not (np.isfinite(x).all() and np.isfinite(z).all()):
            raise ValueError("non-finite values in CSV polyline")
        if (x < -1e-9).any():
            raise ValueError("x must be >= 0 (right-half profile only)")
        if prepend_pillar:
            pillar_x = np.array([0.0, Rxy_p])
            pillar_z = np.array([-Hp, -Hp])
            if not (np.isclose(x[0], Rxy_p) and z[0] <= 0):
                pillar_x = np.concatenate([pillar_x, [Rxy_p]])
                pillar_z = np.concatenate([pillar_z, [z[0]]])
            x = np.concatenate([pillar_x, x])
            z = np.concatenate([pillar_z, z])
        if abs(x[0]) > 1e-6:
            raise ValueError(
                f"first point must be on the axis (x≈0); got x={x[0]:.4f}. "
                "Enable 'prepend pillar' if your CSV only contains the lobe."
            )
        if abs(x[-1]) > 1e-6:
            raise ValueError(
                f"last point must be on the axis (x≈0); got x={x[-1]:.4f}"
            )
        return x, z

    def load_svg_bytes(data, Hp=220.0, units="um", samples=200, path_index=0,
                       flip_y=True, prepend_pillar=False, Rxy_p=25.0):
        from svgpathtools import svgstr2paths
        scale = _UM_PER_UNIT[units]
        text = data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else data
        paths, _attrs = svgstr2paths(text)
        if not paths:
            raise ValueError("no <path> elements found in SVG")
        if path_index >= len(paths):
            raise IndexError(
                f"path_index={path_index} out of range (file has {len(paths)} paths)"
            )
        svg_path = paths[path_index]
        ts = np.linspace(0.0, 1.0, samples)
        pts = np.array([svg_path.point(t) for t in ts])
        x = pts.real.astype(float)
        y = pts.imag.astype(float)
        if flip_y:
            y = -y
        x -= x.min()
        x *= scale
        z = y * scale
        if x[0] > 1e-3:
            i0 = int(np.argmin(x))
            x = np.concatenate([x[i0:], x[:i0]])
            z = np.concatenate([z[i0:], z[:i0]])
        on_axis = np.where(x < 1e-3)[0]
        if len(on_axis) >= 2:
            x = x[on_axis[0]:on_axis[-1] + 1]
            z = z[on_axis[0]:on_axis[-1] + 1]
            x[0] = 0.0
            x[-1] = 0.0
        if prepend_pillar:
            pillar_x = np.array([0.0, Rxy_p, Rxy_p])
            pillar_z = np.array([-Hp, -Hp, z[0]])
            x = np.concatenate([pillar_x, x])
            z = np.concatenate([pillar_z, z])
        if abs(x[0]) > 1e-3 or abs(x[-1]) > 1e-3:
            raise ValueError(
                f"polyline endpoints not on axis (start x={x[0]:.4f}, "
                f"end x={x[-1]:.4f}). SVG path must start and end on the axis."
            )
        return x, z

    def load_csv_outline_bytes(data, units="um"):
        """Full closed 2D outline loader (no axis contract, no pillar)."""
        scale = _UM_PER_UNIT[units]
        text = data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else data
        reader = _csv.reader(_io.StringIO(text))
        rows_raw = [r for r in reader if r and not r[0].startswith("#")]
        if not rows_raw:
            raise ValueError("CSV is empty")
        if _looks_like_header(rows_raw[0]):
            rows_raw = rows_raw[1:]
        x = np.array([float(r[0]) for r in rows_raw]) * scale
        z = np.array([float(r[1]) for r in rows_raw]) * scale
        if not (np.isfinite(x).all() and np.isfinite(z).all()):
            raise ValueError("non-finite values in CSV outline")
        return x, z

    def load_svg_outline_bytes(data, units="um", samples=300, flip_y=True):
        """Full closed 2D outline loader across ALL <path> elements."""
        from svgpathtools import svgstr2paths
        scale = _UM_PER_UNIT[units]
        text = data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else data
        paths, _attrs = svgstr2paths(text)
        if not paths:
            raise ValueError("no <path> elements found in SVG")
        outlines = []
        for svg_path in paths:
            ts = np.linspace(0.0, 1.0, samples)
            pts = np.array([svg_path.point(t) for t in ts])
            x = pts.real.astype(float)
            y = pts.imag.astype(float)
            if flip_y:
                y = -y
            outlines.append((x * scale, y * scale))
        return outlines  # list of (x, z) per path

    return (
        load_csv_bytes, load_csv_outline_bytes,
        load_svg_bytes, load_svg_outline_bytes,
    )


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

    def top_body_vertices(half_x, half_z, P, plate_thickness, initial_gap,
                          bottom_z_apex=None):
        Hp = -float(half_z[0])
        z_apex_top = float(np.asarray(half_z).max())
        if bottom_z_apex is None:
            bottom_z_apex = z_apex_top
        mirror_offset = float(bottom_z_apex) + z_apex_top
        inner = Hp + mirror_offset + initial_gap
        outer = inner + plate_thickness
        top_x_right = P / 2.0 + np.asarray(half_x)
        top_z = mirror_offset + initial_gap - np.asarray(half_z)
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

    def outline_body_layout(outlines_xz, P, plate_thickness, mirror_z=None,
                             plate_above=False):
        """Build a body from one-or-more closed 2D outlines + auto plate.

        - `outlines_xz`: list of (x, z) numpy arrays (closed loops).
        - `P`: cell pitch in µm.
        - `plate_thickness`: thickness of the auto-added plate.
        - `mirror_z`: if not None, z-coordinate to mirror the outline
          across (used to flip the top body so its 'bottom' faces the
          bottom body).
        - `plate_above`: when True, the plate sits ABOVE the body
          (top-body convention).

        Each outline is translated so the population is x-centered at
        x=P/2 (the cell midline). The plate spans the full pitch and
        sits directly below (or above) the outline group.
        """
        if not outlines_xz:
            raise ValueError("outline list is empty")

        # Combine bounds across all outlines for centering.
        all_x = np.concatenate([np.asarray(x) for x, _z in outlines_xz])
        all_z = np.concatenate([np.asarray(z) for _x, z in outlines_xz])
        x_center = 0.5 * (all_x.min() + all_x.max())
        x_shift = P / 2.0 - x_center
        z_min = float(all_z.min())
        z_max = float(all_z.max())

        # Move outlines so their z floor is at 0 (then plate below at z<0).
        z_shift = -z_min
        if mirror_z is not None:
            # Mirror across mirror_z: z' = 2*mirror_z - z. Then shift.
            pass  # handled below per-outline

        body_polygons = []
        for x_arr, z_arr in outlines_xz:
            xs = np.asarray(x_arr, dtype=float) + x_shift
            zs = np.asarray(z_arr, dtype=float) + z_shift
            if mirror_z is not None:
                # Reflect around the supplied z (which is in already-shifted coords)
                zs = 2.0 * mirror_z - zs
                # Reversing keeps CCW orientation after reflection.
                xs = xs[::-1]; zs = zs[::-1]
            body_polygons.append((xs, zs))

        z_body_lo = min(z.min() for _x, z in body_polygons)
        z_body_hi = max(z.max() for _x, z in body_polygons)
        if plate_above:
            plate_inner = z_body_hi
            plate_outer = z_body_hi + plate_thickness
            plate_face_tag = "top_disp"
            body_face_tag = "top_contact"
        else:
            plate_inner = z_body_lo
            plate_outer = z_body_lo - plate_thickness
            plate_face_tag = "bottom_fixed"
            body_face_tag = "bottom_contact"

        # Build a single combined polygon. For multiple disjoint
        # outlines, we connect them through the plate (visually one
        # body; the contact tags isolate the surfaces).
        # Simpler approach for now: emit each outline separately, plus
        # one plate polygon. Caller (gmsh assembly) handles multiple
        # surfaces. The marimo preview just plots them.
        return {
            "plate": [
                (0.0, plate_outer, plate_face_tag),
                (P,   plate_outer, "sym_right"),
                (P,   plate_inner, "sym_right"),
                (0.0, plate_inner, "sym_left"),
            ],
            "outlines": body_polygons,
            "body_face_tag": body_face_tag,
            "z_lo": min(z_body_lo, plate_outer, plate_inner),
            "z_hi": max(z_body_hi, plate_outer, plate_inner),
            "z_apex": z_body_hi,
            "z_floor": z_body_lo,
        }

    return outline_body_layout, bottom_body_vertices, top_body_vertices


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
def _(mo):
    # One global length-unit selector. All geometry sliders, uploaded
    # CSV/SVG values, and the YAML preview are interpreted in this unit
    # — internal pipeline still uses µm.
    units = mo.ui.dropdown(
        options=["um", "mm", "cm", "m"], value="um",
        label="Length units",
    )
    return (units,)


@app.cell
def _(SHAPES, mo):
    # Bottom-body profile source. `builtin` keeps the original dropdown
    # workflow; `csv`/`svg` swap in a file uploader.
    bottom_kind = mo.ui.dropdown(
        options=["builtin", "csv", "svg"], value="builtin",
        label="Bottom profile source",
    )
    shape = mo.ui.dropdown(options=sorted(SHAPES.keys()), value="sphere", label="Bottom shape")
    bottom_file = mo.ui.file(filetypes=[".csv", ".svg"], kind="area", label="Bottom CSV/SVG")
    bottom_prepend = mo.ui.checkbox(value=False, label="Prepend pillar to bottom lobe")
    bottom_flip_z = mo.ui.checkbox(value=False, label="Flip bottom body in z")

    asymmetric = mo.ui.checkbox(value=False, label="Asymmetric pair (different top body)")

    top_kind = mo.ui.dropdown(
        options=["builtin", "csv", "svg"], value="builtin", label="Top profile source",
    )
    top_shape = mo.ui.dropdown(options=sorted(SHAPES.keys()), value="cap", label="Top shape")
    top_file = mo.ui.file(filetypes=[".csv", ".svg"], kind="area", label="Top CSV/SVG")
    top_prepend = mo.ui.checkbox(value=False, label="Prepend pillar to top lobe")
    top_flip_z = mo.ui.checkbox(value=False, label="Flip top body in z")

    direction = mo.ui.radio(
        options=["down (push-in)", "up (retention)"],
        value="down (push-in)",
        label="Direction",
    )
    return (
        asymmetric, bottom_file, bottom_flip_z, bottom_kind, bottom_prepend,
        direction, shape, top_file, top_flip_z, top_kind, top_prepend, top_shape,
    )


@app.cell
def _(RECOMMENDED, mo, shape, units):
    # Scale factor: divide µm values by this to express them in user units.
    _UM_PER_UNIT = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1e6}
    _to_unit = 1.0 / _UM_PER_UNIT[units.value]
    rec_ig, rec_disp, rec_steps = RECOMMENDED[shape.value]
    _u = units.value
    initial_gap = mo.ui.number(
        start=-50000.0 * _to_unit, stop=50000.0 * _to_unit,
        step=max(0.5 * _to_unit, 1e-9),
        value=float(rec_ig) * _to_unit,
        label=f"Initial gap ({_u})  —  negative = bodies overlap at start",
    )
    disp = mo.ui.number(
        start=0.5 * _to_unit, stop=100000.0 * _to_unit,
        step=max(0.5 * _to_unit, 1e-9),
        value=float(rec_disp) * _to_unit,
        label=f"Total displacement ({_u})",
    )
    steps = mo.ui.number(
        start=1, stop=5000, step=5, value=int(rec_steps),
        label="Newton increments (~0.2 µm/step recommended)",
    )
    return disp, initial_gap, steps


@app.cell
def _(
    asymmetric, bottom_file, bottom_kind, bottom_prepend,
    direction, disp, initial_gap, mo, shape, steps, units,
    top_file, top_kind, top_prepend, top_shape,
    bottom_flip_z, top_flip_z,
):
    def _body_widget(kind_picker, shape_picker, file_picker, prepend_picker, flip_picker):
        if kind_picker.value == "builtin":
            return mo.vstack([shape_picker, flip_picker])
        return mo.vstack([file_picker, mo.hstack([prepend_picker, flip_picker])])

    bottom_block = mo.vstack([
        bottom_kind,
        _body_widget(bottom_kind, shape, bottom_file, bottom_prepend, bottom_flip_z),
    ])
    top_block = mo.vstack([
        top_kind,
        _body_widget(top_kind, top_shape, top_file, top_prepend, top_flip_z),
    ])

    mo.vstack([
        mo.hstack([units, direction]),
        asymmetric,
        mo.hstack([bottom_block, top_block]) if asymmetric.value else bottom_block,
        initial_gap,
        disp,
        steps,
    ])
    return


@app.cell
def _(mo, units):
    _UM_PER_UNIT = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1e6}
    _to_unit = 1.0 / _UM_PER_UNIT[units.value]
    _u = units.value
    pitch = mo.ui.number(
        start=10.0 * _to_unit, stop=1_000_000.0 * _to_unit,
        step=max(1.0 * _to_unit, 1e-9), value=198.0 * _to_unit,
        label=f"Pitch ({_u})",
    )
    Hp = mo.ui.number(
        start=0.0, stop=1_000_000.0 * _to_unit,
        step=max(1.0 * _to_unit, 1e-9), value=220.0 * _to_unit,
        label=f"Pillar height Hp ({_u})",
    )
    plate = mo.ui.number(
        start=1.0 * _to_unit, stop=10000.0 * _to_unit,
        step=max(1.0 * _to_unit, 1e-9), value=30.0 * _to_unit,
        label=f"Plate thickness ({_u})",
    )
    mesh_h = mo.ui.number(
        start=0.1 * _to_unit, stop=500.0 * _to_unit,
        step=max(0.5 * _to_unit, 1e-9), value=4.0 * _to_unit,
        label=f"Mesh size ({_u})",
    )
    bottom_x_offset = mo.ui.number(
        start=-1_000_000.0 * _to_unit, stop=1_000_000.0 * _to_unit,
        step=max(0.1 * _to_unit, 1e-9), value=0.0,
        label=f"Bottom body x-offset ({_u}, negative = left)",
    )
    top_x_offset = mo.ui.number(
        start=-1_000_000.0 * _to_unit, stop=1_000_000.0 * _to_unit,
        step=max(0.1 * _to_unit, 1e-9), value=0.0,
        label=f"Top body x-offset ({_u}, negative = left)",
    )
    return Hp, bottom_x_offset, mesh_h, pitch, plate, top_x_offset


@app.cell
def _(Hp, bottom_x_offset, disp, initial_gap, mesh_h, pitch, plate,
       top_x_offset, units):
    """Convert user-unit inputs to µm for the rest of the pipeline."""
    _UM_PER_UNIT = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1e6}
    _f = _UM_PER_UNIT[units.value]
    pitch_um = float(pitch.value) * _f
    Hp_um = float(Hp.value) * _f
    plate_um = float(plate.value) * _f
    initial_gap_um = float(initial_gap.value) * _f
    disp_um = float(disp.value) * _f
    mesh_um = float(mesh_h.value) * _f
    bottom_x_offset_um = float(bottom_x_offset.value) * _f
    top_x_offset_um = float(top_x_offset.value) * _f
    return (
        Hp_um, bottom_x_offset_um, disp_um, initial_gap_um,
        mesh_um, pitch_um, plate_um, top_x_offset_um,
    )


@app.cell
def _(Hp, bottom_x_offset, mesh_h, mo, pitch, plate, top_x_offset):
    mo.accordion({
        "Geometry (advanced)": mo.vstack([
            pitch, Hp, plate, mesh_h,
            mo.hstack([bottom_x_offset, top_x_offset]),
        ]),
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
    Hp_um, SHAPES, asymmetric, bottom_file, bottom_kind, bottom_prepend,
    load_csv_bytes, load_csv_outline_bytes, load_svg_bytes, load_svg_outline_bytes,
    mo, shape, units,
    top_file, top_kind, top_prepend, top_shape,
):
    """Resolve each body's geometry. Tries right-half loaders first;
    if they fail (e.g., uploaded file is a full closed cross-section
    that doesn't touch the axis), falls back to outline-mode loaders
    that accept any closed 2D polygon(s).

    Returns:
      - (half_x, half_z) when mode == "right_half" (legacy axisymmetric)
      - list of (x, z) outlines when mode == "outline"
    """
    profile_error = None

    def _try_right_half(kind, uploaded_file, prepend):
        entry = uploaded_file.value[0]
        contents = entry.contents
        name = entry.name.lower()
        is_svg = name.endswith(".svg") or kind == "svg"
        loader = load_svg_bytes if is_svg else load_csv_bytes
        x, z = loader(contents, Hp=Hp_um, units=units.value,
                      prepend_pillar=prepend.value)
        # If the loader returned a closed loop (start == end as a 2D
        # point), it's actually a full outline that happens to satisfy
        # the axis check trivially. Force the outline fallback so the
        # downstream code treats it as a closed body, not a right-half.
        if abs(x[0] - x[-1]) < 1e-6 and abs(z[0] - z[-1]) < 1e-6:
            raise ValueError("closed loop — use outline mode")
        return x, z

    def _try_outline(kind, uploaded_file):
        entry = uploaded_file.value[0]
        contents = entry.contents
        name = entry.name.lower()
        is_svg = name.endswith(".svg") or kind == "svg"
        if is_svg:
            return load_svg_outline_bytes(contents, units=units.value)
        x, z = load_csv_outline_bytes(contents, units=units.value)
        return [(x, z)]  # CSV → single polygon

    def _resolve(kind, builtin_name, uploaded_file, prepend):
        """Return (mode, data). data shape depends on mode."""
        if kind == "builtin":
            return "right_half", SHAPES[builtin_name](Hp=Hp_um)
        if not uploaded_file.value:
            return None, None
        # Try right-half first (back-compat). Fall back to outline.
        try:
            return "right_half", _try_right_half(kind, uploaded_file, prepend)
        except Exception:
            return "outline", _try_outline(kind, uploaded_file)

    try:
        bot_mode, bot_data = _resolve(bottom_kind.value, shape.value, bottom_file, bottom_prepend)
    except Exception as e:  # noqa: BLE001
        profile_error = f"Bottom profile error: {e}"
        bot_mode, bot_data = "right_half", SHAPES["sphere"](Hp=Hp_um)

    if bot_data is None:
        profile_error = "Upload a bottom CSV/SVG, or switch source back to 'builtin'."
        bot_mode, bot_data = "right_half", SHAPES["sphere"](Hp=Hp_um)

    if asymmetric.value:
        try:
            top_mode, top_data = _resolve(top_kind.value, top_shape.value,
                                          top_file, top_prepend)
        except Exception as e:  # noqa: BLE001
            profile_error = (profile_error + " | " if profile_error else "") + f"Top profile error: {e}"
            top_mode, top_data = "right_half", SHAPES["cap"](Hp=Hp_um)
        if top_data is None:
            profile_error = (profile_error + " | " if profile_error else "") + "Upload a top CSV/SVG."
            top_mode, top_data = "right_half", SHAPES["cap"](Hp=Hp_um)
    else:
        top_mode, top_data = bot_mode, bot_data

    return bot_data, bot_mode, profile_error, top_data, top_mode


@app.cell
def _(mo, profile_error):
    mo.callout(profile_error, kind="warn") if profile_error else None
    return


@app.cell
def _(
    asymmetric, bot_data, bot_mode, bottom_body_vertices,
    bottom_flip_z, bottom_x_offset_um, dir_value, disp_um, initial_gap_um,
    np, pitch_um, plate_um, plt, shape, top_body_vertices, top_data, top_mode,
    top_flip_z, top_shape, top_x_offset_um, units,
):
    """Two-panel preview supporting both right-half and outline modes."""
    _UM_PER_UNIT = {"um": 1.0, "mm": 1_000.0, "cm": 10_000.0, "m": 1e6}
    _to_unit = 1.0 / _UM_PER_UNIT[units.value]
    _u = units.value

    def _plate_polygon(P, z_inner, z_outer):
        """Cell-spanning plate rectangle. CCW when z_outer < z_inner."""
        return (np.array([0.0, P, P, 0.0]),
                np.array([z_outer, z_outer, z_inner, z_inner]))

    def _outlines_centered(outlines, P):
        """X-center a list of (x, z) outlines on x=P/2."""
        all_x = np.concatenate([np.asarray(x) for x, _z in outlines])
        x_shift = P / 2.0 - 0.5 * (float(all_x.min()) + float(all_x.max()))
        return [(np.asarray(x) + x_shift, np.asarray(z)) for x, z in outlines]

    def _flip_z_outlines(outlines):
        """Negate z and reverse vertex order (keeps CCW after flip)."""
        return [(np.asarray(x)[::-1], -np.asarray(z)[::-1]) for x, z in outlines]

    # --- BOTTOM body ---
    if bot_mode == "right_half":
        half_x, half_z = bot_data
        if bottom_flip_z.value:
            # Negate z (apex now negative) and reverse order so half_z[0]
            # is still the "base" point that bottom_body_vertices uses.
            half_x = np.asarray(half_x)[::-1]
            half_z = -np.asarray(half_z)[::-1]
        bvx, bvz, _ = bottom_body_vertices(half_x, half_z, pitch_um, plate_um)
        bottom_polys = [(np.asarray(bvx) + bottom_x_offset_um, np.asarray(bvz))]
        bot_z_apex = float(np.asarray(half_z).max())
    else:
        bot_outlines = _flip_z_outlines(bot_data) if bottom_flip_z.value else bot_data
        centered_b = _outlines_centered(bot_outlines, pitch_um)
        all_zb = np.concatenate([np.asarray(z) for _x, z in centered_b])
        z_shift = -float(all_zb.min())
        positioned_b = [(x + bottom_x_offset_um, z + z_shift) for x, z in centered_b]
        bot_z_apex = float(max(z.max() for _x, z in positioned_b))
        bottom_polys = [_plate_polygon(pitch_um, 0.0, -plate_um)]
        bottom_polys.extend(positioned_b)

    # --- TOP body ---
    if top_mode == "right_half":
        bottom_z_apex_arg = bot_z_apex if asymmetric.value else None
        top_half_x, top_half_z = top_data
        if top_flip_z.value:
            top_half_x = np.asarray(top_half_x)[::-1]
            top_half_z = -np.asarray(top_half_z)[::-1]
        tvx, tvz_init, _ = top_body_vertices(
            top_half_x, top_half_z, pitch_um, plate_um, initial_gap_um,
            bottom_z_apex=bottom_z_apex_arg,
        )
        top_polys_init = [(np.asarray(tvx) + top_x_offset_um, np.asarray(tvz_init))]
    else:
        top_outlines = _flip_z_outlines(top_data) if top_flip_z.value else top_data
        centered_t = _outlines_centered(top_outlines, pitch_um)
        all_zt = np.concatenate([np.asarray(z) for _x, z in centered_t])
        z_max_t = float(all_zt.max())
        # Mirror upside down (so the lobe faces the pin), then translate.
        mirrored = [(np.asarray(x)[::-1], -np.asarray(z)[::-1]) for x, z in centered_t]
        floor_after_mirror = -z_max_t
        dz = (bot_z_apex + initial_gap_um) - floor_after_mirror
        positioned_t = [(x + top_x_offset_um, z + dz) for x, z in mirrored]
        top_z_hi = float(max(z.max() for _x, z in positioned_t))
        top_polys_init = [_plate_polygon(pitch_um, top_z_hi, top_z_hi + plate_um)]
        top_polys_init.extend(positioned_t)

    sign = -1.0 if dir_value == "down" else +1.0
    top_polys_final = [(x, z + sign * disp_um) for x, z in top_polys_init]

    fig, axes = plt.subplots(1, 2, figsize=(10, 8), dpi=110)
    all_x = np.concatenate([x for x, _z in bottom_polys + top_polys_init + top_polys_final])
    all_z = np.concatenate([z for _x, z in bottom_polys + top_polys_init + top_polys_final])
    x_lo = min(all_x.min(), 0.0); x_hi = max(all_x.max(), pitch_um)
    z_lo = float(all_z.min()); z_hi = float(all_z.max())
    x_pad = 0.05 * (x_hi - x_lo + 1.0)
    z_pad = 0.02 * (z_hi - z_lo + 1.0)

    for ax, top_polys, title in [
        (axes[0], top_polys_init,  "initial   (top plate disp = 0)"),
        (axes[1], top_polys_final, f"final     (top plate disp = {(sign * disp_um) * _to_unit:+.3g} {_u})"),
    ]:
        for x, z in bottom_polys:
            ax.fill(x * _to_unit, z * _to_unit,
                    facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.0)
        for x, z in top_polys:
            ax.fill(x * _to_unit, z * _to_unit,
                    facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.0)
        ax.axvline(0.0,                 linestyle="--", color="#888", linewidth=0.8)
        ax.axvline(pitch_um * _to_unit, linestyle="--", color="#888", linewidth=0.8)
        ax.set_aspect("equal")
        ax.set_xlim((x_lo - x_pad) * _to_unit, (x_hi + x_pad) * _to_unit)
        ax.set_ylim((z_lo - z_pad) * _to_unit, (z_hi + z_pad) * _to_unit)
        ax.grid(alpha=0.2)
        ax.set_xlabel(f"x ({_u})")
        ax.set_ylabel(f"z ({_u})")
        ax.set_title(title, fontsize=10)

    bot_label = shape.value if bot_mode == "right_half" else "outline"
    top_label = top_shape.value if top_mode == "right_half" else "outline"
    pair_label = f"{bot_label} ↔ {top_label}" if asymmetric.value else bot_label
    fig.suptitle(
        f"{pair_label} — IG = {initial_gap_um * _to_unit:+g} {_u}, "
        f"DISP = {disp_um * _to_unit:g} {_u}, "
        f"pitch = {pitch_um * _to_unit:g} {_u}",
        fontsize=10,
    )
    fig.tight_layout()
    fig
    return


@app.cell
def _(
    Hp_um, asymmetric, bottom_file, bottom_kind, bottom_prepend,
    contact_mode, dir_value, disp_um, fric, gamma_scale, initial_gap_um, io,
    mesh_um, newton_atol, newton_rtol, out_dir, pitch_um, plate_um, shape,
    steps, top_file, top_kind, top_prepend, top_shape, units, yaml,
):
    """Build the YAML config from current UI state. All geometry values
    are converted to µm (the schema's canonical unit). Uploaded
    profiles keep the user's unit declaration so the local CLI loads
    them with the right scale."""

    def _profile_block(kind, builtin_name, uploaded, prepend):
        if kind == "builtin":
            return {"kind": "builtin", "name": builtin_name, "units": "um"}
        if not uploaded.value:
            return {"kind": kind, "path": "profiles/MISSING_UPLOAD",
                    "units": units.value, "prepend_pillar": bool(prepend.value)}
        fname = uploaded.value[0].name
        return {"kind": kind, "path": f"profiles/{fname}",
                "units": units.value, "prepend_pillar": bool(prepend.value)}

    cfg_dict = {
        "geometry": {
            "pitch_um": float(pitch_um),
            "Hp_um": float(Hp_um),
            "plate_um": float(plate_um),
            "initial_gap_um": float(initial_gap_um),
        },
        "profile": _profile_block(
            bottom_kind.value, shape.value, bottom_file, bottom_prepend,
        ),
    }
    if asymmetric.value:
        cfg_dict["top_profile"] = _profile_block(
            top_kind.value, top_shape.value, top_file, top_prepend,
        )
    cfg_dict.update({
        "mesh": {"characteristic_size_um": float(mesh_um)},
        "solver": {
            "steps": int(steps.value),
            "disp_um": float(disp_um),
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
    })
    _buf = io.StringIO()
    yaml.safe_dump(cfg_dict, _buf, sort_keys=False, default_flow_style=False)
    yaml_text = _buf.getvalue()
    return (yaml_text,)


@app.cell
def _(mo, yaml_text):
    mo.md(f"### YAML preview\n\n```yaml\n{yaml_text}\n```")
    return


@app.cell
def _(
    asymmetric, bottom_file, bottom_kind, dir_value, io, mo, shape,
    top_file, top_kind, top_shape, yaml_text,
):
    """Download cell.

    With no file uploads: a plain YAML download (back-compat).
    With at least one CSV/SVG upload: a ZIP bundle containing the YAML
    plus `profiles/<name>` files, matching the paths emitted in the
    YAML preview. The user unzips the bundle and runs the CLI from
    the unzipped directory.
    """
    import zipfile

    pair = f"{shape.value}_vs_{top_shape.value}" if asymmetric.value else shape.value
    base = f"{pair}_{dir_value}"

    files_to_bundle = []
    if bottom_kind.value != "builtin" and bottom_file.value:
        files_to_bundle.append(bottom_file.value[0])
    if asymmetric.value and top_kind.value != "builtin" and top_file.value:
        files_to_bundle.append(top_file.value[0])

    if files_to_bundle:
        _buf = io.BytesIO()
        with zipfile.ZipFile(_buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(f"{base}.yaml", yaml_text)
            for entry in files_to_bundle:
                zf.writestr(f"profiles/{entry.name}", entry.contents)
        fname = f"{base}.zip"
        download = mo.download(
            data=_buf.getvalue(),
            filename=fname,
            mimetype="application/zip",
            label=f"Download {fname} (bundle: YAML + uploaded profile{'s' if len(files_to_bundle) > 1 else ''})",
        )
    else:
        fname = f"{base}.yaml"
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
