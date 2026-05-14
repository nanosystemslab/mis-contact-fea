"""Streamlit-based 2D alignment GUI.

Replaces the trial-and-error loop of editing `preview_positions.py`,
re-rendering, eyeballing, repeating. Run with:

    streamlit run src/mis_contact_fea/gui/align_2d.py
    # or use the convenience wrapper:
    ./scripts/launch_gui.sh

What you can do:
- Pick a built-in profile (sphere, oblate, ...).
- Adjust pitch, initial gap, displacement, steps, direction with sliders.
- See a live two-panel preview: initial layout (disp=0) and final
  layout (disp=DISP) side by side.
- Tweak Newton / contact-solver knobs in an "advanced" expander.
- Click "Save config" to download a YAML you can immediately feed to
  `./scripts/run_config.sh <path>.yaml`.
- Or load an existing YAML to use as the starting point and edit from
  there.
"""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import yaml

from mis_contact_fea.config import (
    GeometryConfig,
    MeshConfig,
    OutputConfig,
    ProfileConfig,
    SimulationConfig,
    SolverConfig,
)
from mis_contact_fea.mesh.layout import bottom_body_vertices, top_body_vertices
from mis_contact_fea.profiles import SHAPES


# Tuned starting points per shape (matched to examples/01_seven_builtins).
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


def render_preview(shape: str, initial_gap: float, disp: float, pitch: float,
                   plate: float, Hp: float, direction: str) -> plt.Figure:
    """Two-panel matplotlib preview: bodies at start and end of the run."""
    half_x, half_z = SHAPES[shape](Hp=Hp)
    bvx, bvz, _ = bottom_body_vertices(half_x, half_z, pitch, plate)
    tvx, tvz_init, _ = top_body_vertices(half_x, half_z, pitch, plate, initial_gap)
    bvx, bvz = np.asarray(bvx), np.asarray(bvz)
    tvx, tvz_init = np.asarray(tvx), np.asarray(tvz_init)

    # Sign convention matches the solver: direction='down' means top
    # plate moves -z; 'up' means +z.
    sign = -1.0 if direction == "down" else +1.0
    tvz_final = tvz_init + sign * disp

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=120)
    z_lo = min(bvz.min(), tvz_init.min(), tvz_final.min()) - 10
    z_hi = max(bvz.max(), tvz_init.max(), tvz_final.max()) + 10

    for ax, tvz, title in [
        (axes[0], tvz_init,  f"initial   (top plate disp = 0)"),
        (axes[1], tvz_final, f"final     (top plate disp = {sign*disp:+.2f} µm)"),
    ]:
        ax.fill(bvx, bvz, facecolor="#A8C8E8", edgecolor="#1F4F8B", linewidth=1.2)
        ax.fill(tvx, tvz, facecolor="#F5C58C", edgecolor="#9C4A14", linewidth=1.2)
        ax.axvline(0,     linestyle="--", color="#888", linewidth=0.8)
        ax.axvline(pitch, linestyle="--", color="#888", linewidth=0.8)
        ax.set_aspect("equal")
        ax.set_xlim(-15, pitch + 15)
        ax.set_ylim(z_lo, z_hi)
        ax.grid(alpha=0.2)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("z (µm)")
        ax.set_title(title)
    fig.suptitle(f"{shape} — INITIAL_GAP = {initial_gap:+g} µm,  DISP = {disp:g} µm,  pitch = {pitch:g} µm",
                 fontsize=10)
    fig.tight_layout()
    return fig


def build_config(shape: str, geometry_kwargs: dict, mesh_kwargs: dict,
                 solver_kwargs: dict, out_dir: str) -> SimulationConfig:
    return SimulationConfig(
        geometry=GeometryConfig(**geometry_kwargs),
        profile=ProfileConfig(kind="builtin", name=shape),
        mesh=MeshConfig(**mesh_kwargs),
        solver=SolverConfig(**solver_kwargs),
        output=OutputConfig(out_dir=Path(out_dir)),
    )


def main() -> None:
    st.set_page_config(
        page_title="mis-contact-fea — alignment GUI",
        layout="wide",
    )
    st.title("mis-contact-fea — 2D alignment")
    st.caption(
        "Set up `INITIAL_GAP`, `DISP`, and solver parameters visually, "
        "then export a YAML config you can run with "
        "`./scripts/run_config.sh <yaml>`."
    )

    # ---- Sidebar: load existing config (optional starting point) ----
    with st.sidebar:
        st.header("Profile & motion")
        loaded: SimulationConfig | None = None
        uploaded = st.file_uploader("Load existing config (optional)", type=["yaml", "yml"])
        if uploaded is not None:
            try:
                loaded = SimulationConfig.model_validate(yaml.safe_load(uploaded))
                st.success(f"Loaded {uploaded.name}")
            except Exception as e:  # noqa: BLE001
                st.error(f"Couldn't parse config: {e}")

        # Defaults — prefer values from a loaded config, fall back to
        # the RECOMMENDED table for the chosen shape.
        default_shape = (loaded.profile.name if loaded and loaded.profile.kind == "builtin"
                         else "sphere")
        shape = st.selectbox("Shape", sorted(SHAPES.keys()),
                             index=sorted(SHAPES.keys()).index(default_shape))

        rec_ig, rec_disp, rec_steps = RECOMMENDED[shape]
        ig_default   = loaded.geometry.initial_gap_um if loaded else rec_ig
        disp_default = loaded.solver.disp_um         if loaded else rec_disp
        steps_default = loaded.solver.steps          if loaded else rec_steps
        dir_default  = (loaded.solver.direction if loaded else "down")
        pitch_default = loaded.geometry.pitch_um if loaded else 198.0
        Hp_default    = loaded.geometry.Hp_um    if loaded else 220.0
        plate_default = loaded.geometry.plate_um if loaded else 30.0
        mesh_default  = (loaded.mesh.characteristic_size_um if loaded else 4.0)

        direction = st.radio("Direction", ["down (push-in)", "up (retention)"],
                             index=0 if dir_default == "down" else 1,
                             horizontal=True)
        direction_value = "down" if direction.startswith("down") else "up"

        initial_gap = st.slider(
            "Initial gap (µm)  —  negative = bodies overlap at start",
            min_value=-250.0, max_value=50.0, value=float(ig_default), step=0.5,
        )
        disp = st.slider(
            "Total displacement (µm)",
            min_value=0.5, max_value=300.0, value=float(disp_default), step=0.5,
        )
        steps = st.number_input(
            "Newton increments (~0.2 µm/step recommended)",
            min_value=1, max_value=5000, value=int(steps_default), step=5,
        )

        st.markdown("---")
        st.subheader("Geometry (advanced)")
        pitch = st.slider("Pitch (µm)", 50.0, 400.0, float(pitch_default), step=1.0)
        Hp = st.slider("Pillar height Hp (µm)", 80.0, 400.0, float(Hp_default), step=5.0)
        plate = st.slider("Plate thickness (µm)", 5.0, 100.0, float(plate_default), step=1.0)
        mesh_h = st.slider("Mesh size (µm)", 1.0, 20.0, float(mesh_default), step=0.5)

        st.markdown("---")
        st.subheader("Solver (advanced)")
        with st.expander("Newton / contact knobs", expanded=False):
            contact_mode = st.selectbox(
                "Contact mode", ["raytracing", "closest"],
                index=0 if (loaded.solver.contact_mode if loaded else "raytracing") == "raytracing" else 1,
                help="raytracing is more robust for push; closest is preferred for retention.",
            )
            fric = st.number_input(
                "Friction coefficient", 0.0, 2.0,
                float(loaded.solver.fric if loaded else 0.0), step=0.05,
            )
            newton_rtol = st.number_input(
                "Newton rtol", 1e-6, 1.0,
                float(loaded.solver.newton_rtol if loaded else 0.2),
                format="%.4g",
            )
            newton_atol = st.number_input(
                "Newton atol", 1e-6, 1e6,
                float(loaded.solver.newton_atol if loaded else 1000.0),
                format="%.4g",
            )
            gamma_scale = st.number_input(
                "Gamma scale (Nitsche penalty)", 0.001, 100.0,
                float(loaded.solver.gamma_scale if loaded else 0.2),
                format="%.4g",
            )

        st.markdown("---")
        st.subheader("Output")
        default_out_dir = (
            str(loaded.output.out_dir) if loaded else
            f"runs_all_shapes/{shape}" if direction_value == "down" else
            f"runs_retention/{shape}"
        )
        out_dir = st.text_input("Output directory (relative to project root)",
                                value=default_out_dir)

    # ---- Main area: preview + export ----
    col_preview, col_save = st.columns([3, 1])

    with col_preview:
        st.subheader("Layout preview")
        fig = render_preview(shape, initial_gap, disp, pitch, plate, Hp, direction_value)
        st.pyplot(fig, clear_figure=True)

    with col_save:
        st.subheader("Export")
        cfg = build_config(
            shape,
            geometry_kwargs=dict(
                pitch_um=pitch, Hp_um=Hp, plate_um=plate,
                initial_gap_um=initial_gap,
            ),
            mesh_kwargs=dict(characteristic_size_um=mesh_h),
            solver_kwargs=dict(
                steps=int(steps), disp_um=disp, direction=direction_value,
                contact_mode=contact_mode, fric=fric,
                newton_rtol=newton_rtol, newton_atol=newton_atol,
                gamma_scale=gamma_scale,
            ),
            out_dir=out_dir,
        )

        # Render the YAML for both preview and download.
        buf = io.StringIO()
        yaml.safe_dump(cfg.model_dump(mode="json", exclude_none=True),
                       buf, sort_keys=False, default_flow_style=False)
        yaml_text = buf.getvalue()

        default_fname = f"{shape}_{direction_value}.yaml"
        fname = st.text_input("Filename", value=default_fname, key="fname_input")
        st.download_button(
            "Download config",
            data=yaml_text,
            file_name=fname,
            mime="text/yaml",
        )

        st.markdown("**Or write to repo:**")
        save_dir = st.text_input(
            "Save path (relative to project root)",
            value=f"examples/01_seven_builtins/{fname}",
            key="save_path_input",
        )
        if st.button("Save to file"):
            try:
                project_root = Path(__file__).resolve().parents[3]
                target = project_root / save_dir
                target.parent.mkdir(parents=True, exist_ok=True)
                cfg.to_yaml(target)
                st.success(f"Wrote {target.relative_to(project_root)}")
            except Exception as e:  # noqa: BLE001
                st.error(f"Couldn't save: {e}")

    st.markdown("---")
    st.subheader("YAML preview")
    st.code(yaml_text, language="yaml")
    st.caption(
        "Run with:  `./scripts/run_config.sh "
        f"{save_dir if 'save_dir' in dir() else default_fname}`"
    )


if __name__ == "__main__":
    main()
