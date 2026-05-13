# mis-contact-fea

[![status](https://img.shields.io/badge/status-alpha-orange)](https://github.com/nanosystemslab/mis-contact-fea)
[![python](https://img.shields.io/badge/python-3.10+-blue)](https://www.python.org/downloads/)
[![license](https://img.shields.io/badge/license-GPL--3.0-blue)](LICENSE)
[![DOLFINx](https://img.shields.io/badge/DOLFINx-0.9.0-blue)](https://github.com/FEniCS/dolfinx)
[![dolfinx_contact](https://img.shields.io/badge/dolfinx__contact-Nitsche-blue)](https://github.com/Wells-Group/asimov-contact)

[![Physics](https://img.shields.io/badge/physics-plane--strain%20FEA%20%2B%20contact-blueviolet)](#how-it-works)
[![Shapes](https://img.shields.io/badge/lobe%20shapes-7-green)](#geometry)
[![Mesh](https://img.shields.io/badge/mesh-Gmsh%20parametric-yellow)](src/build_2d_slice_mesh.py)
[![HPC](https://img.shields.io/badge/HPC-SLURM%20%2B%20Singularity-lightgrey)](hpc/)
[![Container](https://img.shields.io/badge/docker-mattnakamura%2Fdolfinx--contact-2496ED?logo=docker)](https://hub.docker.com/r/mattnakamura/dolfinx-contact)

**Mechanical Interlocking Structure — Contact FEA.** A 2D plane-strain
finite-element pipeline for characterizing the **insertion** vs
**retention** force asymmetry of axisymmetric nanostructure shapes
arranged as a mechanical interlock. Built on
[FEniCSx / dolfinx](https://github.com/FEniCS/dolfinx) +
[dolfinx_contact](https://github.com/Wells-Group/asimov-contact) (Nitsche
frictionless contact), with Gmsh for meshing and PyVista for animation.

The pipeline simulates seven candidate lobe geometries (sphere, oblate,
prolate, cone, cap, capcone, doublecone), each laid out in a periodic
unit cell with pitch P=198 µm and offset P/2, then computes the
displacement-controlled force curves for both **push-in** and
**pull-out** directions. From these it extracts a single figure of
merit — the **retention/insertion ratio** — that ranks each shape by
how strong an asymmetric grip it provides.

The simulation runs as a 2D plane-strain slice through the symmetry
plane of the unit cell, which makes it tractable on a laptop while
remaining quantitatively predictive for periodic patches.

## Repository layout

```
mis-contact-fea/
├── README.md
├── .gitignore
├── .python-version
├── requirements.txt          host-side Python deps (matplotlib, gmsh, pyvista, …)
├── src/                      mesh / solver / analysis Python
│   ├── profiles.py                shape generators (right-half polylines)
│   ├── build_2d_slice_mesh.py     Gmsh → dolfinx mesh
│   ├── run_2d_slice_sim.py        Nitsche contact solver (runs in Docker)
│   ├── animate_pv.py              PyVista GIF + frame export
│   ├── analyze_forces.py          smoothing, chatter detection, R/I ratio
│   ├── compare_forces.py          force-displacement overlays
│   ├── compare_shapes.py          per-shape result-grid composer
│   ├── compute_initial_gaps.py    reference-shape-based gap calibration
│   ├── preview_mesh.py            mesh + BC visualization
│   └── preview_positions.py       initial vs final layout (no simulation)
├── scripts/                  bash drivers for sweeps
│   ├── run_2d_pipeline.sh         single-shape: mesh + solve + animate
│   ├── run_all_shapes.sh          inner loop, called by sweeps below
│   ├── run_final_sweep.sh         push (insertion) sweep, all 7 shapes
│   └── run_retention_sweep.sh     retention sweep (auto-starts from push end-state)
├── docker/                   container build for the contact solver
│   ├── Dockerfile.asimov          Docker image (used on laptop)
│   ├── build_asimov_image.sh
│   ├── dolfinx_v0.9.0_contact.def Apptainer/Singularity .def (used on HPC)
│   └── build_contact_container.sh
├── hpc/                      SLURM scripts for cluster runs
│   ├── build_contact_container.slurm
│   └── run_contact_sim.slurm
├── notebooks/                exploratory Jupyter work
│   └── family_2d_geometries_axisym.ipynb
└── archive/                  legacy 3D pipeline + vendored library
    ├── asimov-contact/            upstream dolfinx_contact source (reference)
    ├── cad_models/                STEP files from the original 3D study
    ├── legacy_3d/                 old 3D Python scripts
    ├── legacy_scripts/            superseded shell drivers
    └── legacy_docs/               old RUNBOOK + usage notes
```

Generated outputs (`runs_all_shapes/`, `runs_retention/`, all top-level
PNGs/CSVs/GIFs) are `.gitignore`'d.

## Quickstart

### 1. Get the Docker image (once)

The solver runs inside a Docker container so dolfinx + dolfinx_contact +
petsc4py + mpi4py don't need to be installed on the host.

**Pull the prebuilt image (recommended, ~30 s):**

```bash
docker pull mattnakamura/dolfinx-contact:v0.9.0
```

Image: [`mattnakamura/dolfinx-contact`](https://hub.docker.com/r/mattnakamura/dolfinx-contact) (~5 GB; dolfinx v0.9.0 + dolfinx_contact built from [asimov-contact](https://github.com/Wells-Group/asimov-contact)).

**Or build from source** (~5 min, if you need to rebuild dolfinx_contact):

```bash
./docker/build_asimov_image.sh
```

This tags the image as both `dolfinx-contact:local` and
`mattnakamura/dolfinx-contact:{v0.9.0,latest}` so a subsequent
`docker push` publishes it.

### 2. Install the host-side Python package

The mesher, animator, and analysis tools run on the host:

```bash
pip install -e .
```

This installs `mis_contact_fea` as an editable package; all the
sweep scripts invoke it as `python -m mis_contact_fea.<module>`.

### 3. Preview layout (no simulation, ~5 s)

```bash
python3 src/preview_positions.py
open preview_positions.png
```

Edit `SHAPE_CONFIGS` at the top of that file (or pass
`--initial-gap-sphere -55 --disp-sphere 85` style overrides) to dial in
each shape's starting gap and total displacement before committing to a
long sweep.

### 4. Run the insertion sweep

```bash
./scripts/run_final_sweep.sh
```

Per-shape parameters live in tables at the top of the script. Override
any of them via env vars, e.g.:

```bash
SHAPES="cone doublecone" ./scripts/run_final_sweep.sh
DISP_sphere=90 STEPS_sphere=450 ./scripts/run_final_sweep.sh
```

Output: `runs_all_shapes/<shape>/{mesh,results,run.log,mesh_preview.png,contact_animation_pv.gif}`.

### 5. Run the retention sweep

```bash
./scripts/run_retention_sweep.sh
```

This automatically starts each shape in the final configuration of its
push-in (i.e. `INITIAL_GAP = push_INITIAL_GAP − full_push_DISP`) and
moves the top plate upward to measure pull-out force. Output goes to
`runs_retention/`.

### 6. Analyze + plot

```bash
python3 src/analyze_forces.py
```

This applies the smoothing pipeline to every force curve, computes
peak insertion and retention forces, extracts the R/I ratio, and writes:

- `force_smoothed_overlay.png` — per-shape raw + smoothed for push & retention
- `force_smoothed_push.png` — all-shape smoothed insertion comparison
- `force_smoothed_retention.png` — all-shape smoothed retention comparison
- `force_ratio.png` — R/I bar chart, sorted high-to-low
- `force_peaks.csv` — summary table

## How it works

### Geometry

A unit cell contains two bodies — bottom (rigid plate + pillar + lobe)
and top (mirror image, translated by P/2). The 2D slice is taken
through the symmetry plane, so symmetry BCs on the left/right cell
edges (`u_x = 0`) are equivalent to spatial periodicity for this offset.

- **Pitch:** 198 µm
- **Offset:** 99 µm (P/2 — true interlock arrangement)
- **Pillar height:** 220 µm
- **Plate thickness:** 30 µm
- **Effective lobe diameter:** 100 µm for all seven shapes (within 4 nm)

The cone and doublecone construction radii are *shrunk* (53.033 µm and
50.849 µm respectively) so that the 2-µm corner fillets don't push
their effective outer diameter past 100 µm — see the calibration math
in `src/profiles.py`.

### Mesh

`src/build_2d_slice_mesh.py` builds a Gmsh 2D triangulated mesh
(~5k–8k nodes per shape, `MESH_SIZE = 4 µm`). Each body is a single
closed polyline (pillar ∪ lobe) discretized as one connected
triangular region. Edges are tagged for Dirichlet BCs (`top_plate`,
`bottom_plate`), contact candidates (`top_lobe`, `bottom_lobe`), and
symmetry (`left_edge`, `right_edge`).

### Solver

`src/run_2d_slice_sim.py` runs inside `dolfinx-contact:local` and
solves linear elasticity (steel-like `E = 200 GPa, ν = 0.3`) under
incremental displacement control on the top plate. Contact is enforced
via Nitsche's method (`dolfinx_contact`), with `--contact-mode
raytracing` for insertion and `closest` for retention.

The solver tolerances are deliberately loose:

```
NEWTON_ATOL=1000
NEWTON_RTOL=0.2
GAMMA_SCALE=0.2
```

The strict defaults (`atol=1e-3, rtol=1e-5, gamma=1`) cause Newton to
throw a `did not converge` error once the bodies separate at the end
of a retention pull-out, even though the residual is at the noise
floor. The loose tolerances let the run complete past disengagement.

### Sweep parameters

| shape       | push IG | push DISP | push STEPS | retention DISP | retention STEPS |
|-------------|---------|-----------|------------|----------------|-----------------|
| sphere      | −55     | 85        | 425        | 50             | 250             |
| oblate      | −30     | 40        | 200        | 25             | 125             |
| prolate     | −80     | 87        | 435        | 55             | 275             |
| cone        | −75     | 30.9      | 155        | 18             | 90              |
| cap         | −65     | 40        | 200        | 40             | 200             |
| capcone     | −110    | 45        | 225        | 45             | 225             |
| doublecone  | −80     | 40        | 200        | 28             | 140             |

All increments are ~0.2 µm/step. Retention displacements are shorter
than push displacements because the bodies fully disengage well before
the push displacement would predict — any further steps just burn
compute on the numerical noise floor.

Prolate push was trimmed to 87 µm to stop just before the deep-engagement
Newton-chatter region (which kicks in around disp ≈ 88 µm).

### Smoothing & chatter handling

Real contact problems can produce two distinct kinds of noise in the
force curve:

1. **Single-point spikes** — the contact integrator emits an outlier
   near a sharp geometric feature or at a topology change.
2. **Multi-point chatter** — Newton oscillates between contact
   configurations in deep engagement, producing sign-flipping residuals
   that look like rapid oscillation in the force.

`src/analyze_forces.py` applies a robust three-stage cleanup:

1. **9-point rolling median** — kills sign-flip chatter (the median of
   `+5000 / −4000 / +3000 / …` is near the underlying trend, not the
   noise).
2. **7-point Hampel filter** with a 2.5×MAD threshold — catches
   residual isolated outliers.
3. **5-point moving average** — cosmetic smoothing for plotting;
   preserves peak shape.

A **chatter mask** is computed from the raw signed force, flagging any
region with ≥4 sign-flips inside a 5-point window. These regions are
shaded red on the per-shape overlay plot, and if the smoothed peak
falls inside a chatter region it is marked `CHATTER` in the summary so
that the value (rolling-median-estimated) is clearly distinguished
from a directly observed peak.

## Output

All forces are reported as **|F_z| in µN per µm of out-of-plane depth**.
Plane-strain assumes an infinitely deep slab; multiply by patch depth
to get a real-world force.

### Last full-sweep results

| rank | shape       | push peak | retention peak | **R/I ratio** |
|------|-------------|----------:|---------------:|--------------:|
| 1    | **cap**     | 578       | 3 182          | **5.50**      |
| 2    | **capcone** | 2 114     | 11 130         | **5.27**      |
| 3    | **oblate**  | 123       | 263            | **2.14**      |
| 4    | prolate     | 16        | 19             | 1.21          |
| 5    | sphere      | 57        | 63             | 1.11          |
| 6    | doublecone  | 45        | 17             | 0.37          |
| 7    | cone        | 6 756     | 1 127          | 0.17          |

R/I > 1: harder to remove than to insert (good interlock).
R/I < 1: easier to remove than to insert (anti-interlock).

The clear winners are **cap** and **capcone** (~5× retention bias),
both of which have a lateral shoulder/undercut that the mating lobe
catches against during pull-out. Pointy taper shapes (cone,
doublecone) are easier to remove than to insert.

## HPC / cluster runs

For longer sweeps or 3D extensions, the `hpc/` directory contains SLURM
submission scripts. They build an Apptainer/Singularity image from
`docker/dolfinx_v0.9.0_contact.def` and run the same Python in batch
mode.

```bash
sbatch hpc/build_contact_container.slurm   # one-time, builds .sif
sbatch hpc/run_contact_sim.slurm           # one sweep job
```

## Caveats & known limits

- **Linear elasticity** only. Plastic deformation is not modeled —
  peak stresses near sharp corners will exceed any realistic material
  yield. The cone/capcone retention numbers in particular would soften
  considerably under an elastic-plastic constitutive law.
- **Frictionless contact.** Friction would raise both insertion and
  retention forces and could shift the R/I ratio.
- **2D plane-strain** is exact for an infinitely-deep periodic slab.
  Real patches have finite depth and edge effects not captured here.
- **Newton chatter** in very deep contact was managed by trimming the
  push displacement; some genuine high-engagement signal may be lost.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE) for the full text.
