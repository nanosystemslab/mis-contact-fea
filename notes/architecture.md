# Architecture — how the pieces fit

## The pipeline (3 stages, two of which run on the host)

```
   YAML config (pydantic-validated, schema in config.py)
        ↓
   ┌─────────────────────┐
   │  1. MESHER  (host)  │   gmsh, meshio
   │  slice_2d.py        │   writes mesh.msh + mesh_tags.json
   └─────────────────────┘
        ↓
   ┌─────────────────────────────────────┐
   │  2. SOLVER  (Docker container)      │   dolfinx + dolfinx_contact
   │  contact_2d.py                      │   writes force_displacement.csv +
   │  run via mattnakamura/dolfinx-      │   contact_results.xdmf
   │  contact:v0.9.0                     │
   └─────────────────────────────────────┘
        ↓
   ┌─────────────────────────────────────┐
   │  3. POSTPROC  (host)                │   pyvista, matplotlib
   │  animate.py + analyze.py            │   writes contact_animation_pv.gif
   └─────────────────────────────────────┘
```

`scripts/run_config.sh` orchestrates all three — reads the YAML on the
host, expands it into shell-safe variables, then invokes each stage
with the right CLI flags.

## Package layout

```
src/mis_contact_fea/
├── config.py             # Pydantic schema for YAML configs
├── profiles/             # Right-half polyline generators / loaders
│   ├── builtins.py       # SHAPES dict: sphere, cone, cap, capcone,
│   │                     # doublecone, oblate, prolate
│   ├── from_csv.py       # CSV → (x_um, z_um) — pure stdlib
│   ├── from_svg.py       # SVG → (x_um, z_um) — needs svgpathtools
│   └── from_step.py      # STEP → (x_um, z_um) — needs gmsh
├── mesh/
│   ├── layout.py         # Pure-numpy polyline-layout helpers
│   │                     # bottom_body_vertices, top_body_vertices,
│   │                     # signed_area, dedup_consecutive
│   └── slice_2d.py       # gmsh assembly + meshio export (host)
├── solver/
│   ├── contact_2d.py     # dolfinx + dolfinx_contact (runs in container)
│   └── contact_3d.py     # 3D version (less mature)
├── postproc/
│   ├── animate.py        # pyvista GIF
│   ├── analyze.py        # force-displacement, energy
│   ├── compare_forces.py # multi-shape comparison
│   ├── preview_mesh.py   # matplotlib mesh preview
│   └── preview_positions.py  # body-only preview before meshing
└── gui/
    ├── align_marimo.py   # Marimo notebook (the WASM-hosted GUI)
    ├── align_2d.py       # Streamlit fallback (legacy)
    └── align_3d.py       # 3D Streamlit GUI
```

## Why the marimo notebook inlines so much

`gui/align_marimo.py` ships duplicated copies of:
- The `SHAPES` dict
- `bottom_body_vertices` / `top_body_vertices`
- CSV / SVG byte-loaders
- The new outline-mode helpers

This is intentional — the notebook also gets exported as a WASM bundle
by `.github/workflows/deploy-gui.yml` and served from GitHub Pages.
Inside Pyodide the package isn't installable, so all the code has to
live in the notebook itself.

`tests/test_marimo_inline_sync.py` extracts the inline functions via
AST walking and compares them against the package versions to catch
drift. If you change the inline copy, the package version must match
(or the test fails).

## Config / schema flow

```
   user-supplied YAML
        ↓
   SimulationConfig.from_yaml(path)  # pydantic validation
        ↓
   ┌──────────────────────────────────────────┐
   │ cfg.geometry   pitch_um, Hp_um, plate_um │
   │                initial_gap_um            │
   │ cfg.profile    ProfileConfig (bottom)    │
   │ cfg.top_profile  ProfileConfig (top)     │
   │ cfg.mesh       characteristic_size_um    │
   │ cfg.solver     steps, disp_um, ...       │
   │ cfg.output     out_dir                   │
   └──────────────────────────────────────────┘
```

The schema treats everything in **µm internally**. The user-unit
display in the Marimo GUI is a presentation layer only.

## Right-half polyline contract

`bottom_body_vertices(half_x, half_z, P, plate_thickness)` and
`top_body_vertices(...)` expect the polyline to obey:

- `half_x[0] = 0` (start on axis)
- `half_z[0] = -Hp` (pillar base — sets `Hp = -half_z[0]` implicitly)
- `half_x[-1] = 0` (end on axis at the apex)
- `half_z[-1] = z_apex` (apex on axis)
- Intermediate points trace the body's right-side surface from
  pillar-shoulder up to apex
- For builtin shapes there's a 2-point "pillar shoulder" transition:
  `half_x[0..1] = [0, Rxy_p]`, both at `half_z = -Hp`

The mesher MIRRORS the right-half across `x=P/2` and ALSO uses the
fact that the body straddles the periodic seam `x=0 ↔ x=P` to create
the full body.

**For our STEP-extracted right-half profiles, `half_z[0] = 0` (not
`-Hp`)** — the pin/socket bodies don't have a separate pillar; their
base sits directly on the plate. This works because `Hp = -half_z[0]
= 0` and the layout helpers handle that case correctly.

## Outline mode (in GUI only, not in CLI yet)

Outline mode treats the body as a **closed 2D polygon** with no
axis-mirror semantics. Used in the Marimo preview for full
cross-sections from STEP files. The GUI's `_outlines_centered()` and
inline outline-body construction code positions the polygon centered
at `x=P/2` with an auto-plate added below (for bottom) or above (for
top).

**Not yet plumbed into the CLI mesher / solver schema.** That's the
biggest pending work item — see `pending_work.md`.

## Length-unit handling

```
   Marimo GUI sliders (user units, e.g. mm)
        ↓ multiply by UM_PER_UNIT[unit]
   internal computation (always µm)
        ↓
   YAML config (always *_um suffix)
        ↓
   CLI pipeline (always µm)
        ↓
   Solver assumes µm / MPa / µN units
```

The hosted Marimo GUI's unit-conversion cell is the single place
where user-units → µm conversion happens. Once converted to `pitch_um`,
`Hp_um`, etc., the rest of the pipeline never sees the original unit.

## Test coverage

```
tests/
├── test_asymmetric_layout.py    # top_body_vertices asymmetric mode
├── test_marimo_inline_sync.py   # GUI inlines vs package versions
├── test_profile_loaders.py      # CSV roundtrip, load_profile dispatch
└── test_profiles_smoke.py       # All 7 builtins generate valid polylines
```

All 20 tests pass as of the snapshot.

## Docker container

Image: `mattnakamura/dolfinx-contact:v0.9.0` (also tagged `:latest`).
Definition: `docker/dolfinx_v0.9.0_contact.def` (Singularity build
file — also produces a `.sif` for HPC use).

The container needs to mount the project dir at `/work` so the solver
can read the mesh and write results. `run_config.sh` does this with
`-v "$PWD:/work" -w /work`.

`PYTHONPATH` inside the container is set so that `mis_contact_fea`
imports work without `pip install` — the container only sees the
`src/` tree via the bind mount.
