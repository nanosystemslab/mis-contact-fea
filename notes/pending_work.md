# Pending work — what's not done yet

Sorted by impact × effort.

## High-impact

### 1. Outline-mode mesher (the big one)

**What's missing**: the CLI mesher (`mesh/slice_2d.py` + `mesh/layout.py`)
only handles right-half axisymmetric profiles. The Marimo GUI handles
full-outline closed 2D bodies for *preview only*. The downloaded YAML
+ CLI pipeline can't run such bodies.

**What's needed**:
1. New `body_vertices_outline(closed_xz, P, plate_thickness)` in
   `mesh/layout.py` that takes a closed 2D polygon and emits gmsh
   vertices + edge tags (with auto-plate underneath).
2. Same for the top body (`body_vertices_outline_top(...)` with
   flipping + positioning relative to `bot_z_apex + initial_gap`).
3. Extend `ProfileConfig` with `mode: Literal["right_half", "outline"]
   = "right_half"`.
4. Extend `slice_2d.build_mesh()` to dispatch on mode.
5. Update `gui/align_marimo.py` YAML emitter to set `mode` correctly.
6. Update `scripts/step_to_profile.py` to default the CSV header
   appropriately (or add a `# mode: outline` comment line).

**Estimated effort**: ~1 day. The trickiest piece is handling
disjoint outlines (the socket has 2 lobes per body — multi-body
unit cell?) — likely punt to "1 polygon per body" for v1.

**Why it matters**: enables non-axisymmetric parts in the simulation
pipeline. Right-half mode is fine for the user's pin/socket because
both are revolves, but future parts won't always be.

### 2. Solver convergence robustness

The pin/socket sim fails at step 6 with default Newton params. We
worked around it (see `convergence_tuning.md`) but real robustness
would come from:

- **Auto-tune Newton tolerances** based on problem scale (compute
  initial residual, set `atol = initial_r * 1e-6` or similar)
- **Adaptive step size** in the solver (if a step takes >K Newton
  iters, bisect; if it converges in <M, double)
- **Expose `newton-max-it`, `relax`** in the schema

**Effort**: ~half day for the schema/CLI plumbing, ~2 days for real
adaptive stepping (would require touching the solver loop).

## Medium-impact

### 3. GUI YAML emitter — fix kind/path mismatch

When the user picks "csv" in the source dropdown but uploads a
`.svg` file (or vice versa), the YAML emits the dropdown value not
the file extension. Result: `kind: csv` with `path: foo.svg` — fails
to load in the CLI.

**Fix**: in the YAML emitter in `align_marimo.py`, derive `kind`
from the uploaded file's extension, not from the dropdown:

```python
def _kind_from(uploaded, kind_picker):
    if not uploaded.value:
        return kind_picker.value
    name = uploaded.value[0].name.lower()
    if name.endswith(".svg"):
        return "svg"
    if name.endswith(".csv"):
        return "csv"
    return kind_picker.value
```

**Effort**: 15 minutes.

### 4. Scale-aware slider defaults when files upload

Right now the geometry sliders keep their nano-scale defaults
(pitch=198 µm, etc.) even when the user uploads mm-scale files.
Add an "infer defaults from upload" path: detect the uploaded
polyline's bounding box and propose `pitch ≈ 1.2 × body_width`,
`disp ≈ z_apex_top × 1.2`, etc.

**Effort**: ~1 hour. Small UX win that prevents the most common
"why does my body look tiny / huge" support question.

### 5. MPI plumbing in scripts/run_config.sh

The solver is MPI-ready (uses `MPI.COMM_WORLD`) but the wrapper
script always runs single-rank. Adding `MPI_RANKS=4 ./scripts/
run_config.sh ...` would 2-3× speed up local runs.

**Edit needed** at line ~86 of `scripts/run_config.sh`:

```bash
MPI_RANKS="${MPI_RANKS:-1}"
docker run --rm -v "$PWD:/work" -w /work -e PYTHONPATH="$CONTAINER_PYTHONPATH" "$IMAGE" \
    bash -c "mpirun -n $MPI_RANKS --allow-run-as-root \
        python3 -m mis_contact_fea.solver.contact_2d ..."
```

**Effort**: 10 minutes. Tested via env var so default behavior
unchanged.

### 6. STEP loader in the package (handle axis-along-x case)

`mis_contact_fea/profiles/from_step.py` exists but assumes the
body's axis-of-revolution is z. For the user's STEP files (axis is
x), it produces garbage. The good logic is now in
`scripts/step_to_profile.py`; should be promoted into the package
loader so STEP files work via `kind: step` in YAML configs.

**Effort**: ~2 hours of refactoring + tests.

## Low-impact / cleanup

### 7. Update Node.js actions to v24 (GitHub deprecation)

Pages deploy workflow warns about Node.js 20 → 24 transition.
Edit `.github/workflows/deploy-gui.yml`:

- `actions/checkout@v4` → `@v5` (or set `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24=true`)
- `actions/setup-python@v5` → `@v6` similarly
- `actions/upload-pages-artifact@v3` and `actions/deploy-pages@v4`

**Effort**: 5 minutes. **Deadline**: June 2nd 2026.

### 8. Marimo notebook formatting warnings

Two pre-existing warnings from `marimo check`:
- Docstring at line 35 expected before marimo import
- Markdown indentation

Cosmetic. Could fix in 5 minutes.

### 9. Pin / socket retention example

`pin_socket_retention.yaml` exists but hasn't been tested end-to-end.
Once the push config converges, run retention and verify the
direction=up / contact-mode=closest path works.

**Effort**: depends on whether retention has its own convergence
issues. ~1 hour to try.

### 10. Outline-mode visualization in the GUI

When uploading outline-mode files, the title still says
`outline ↔ outline` instead of the uploaded filenames. Minor UX —
use `bottom_file.value[0].name` instead.

**Effort**: 10 minutes.

## Known design limitations (not bugs — design decisions to revisit)

- **Plate sits across full cell width** even when the body is
  x-offset. This is consistent with the periodic boundary model
  but means an off-center body looks weird (body sticks out past
  the plate). Could either clip the body or extend the plate.

- **Single body per side per cell**. For socket-with-two-lobes the
  current code picks the largest face. True multi-body simulation
  is the outline-mode work.

- **No friction yet tested at scale**. `fric: 0.0` everywhere.
  Should sweep friction once a converging baseline exists.

- **No 3D outline-mode mesher**. `solid_3d.py` is right-half-only.
  Probably fine — most 2D contact studies are sufficient.
