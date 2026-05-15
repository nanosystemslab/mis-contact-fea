# Known bugs — found but not fixed

## GUI

### YAML emitter: `kind` from dropdown not file extension

`gui/align_marimo.py`, in the YAML-emitter cell, builds `profile.kind`
from `bottom_kind.value` (the dropdown). If the user picks "csv" in
the dropdown but uploads a `.svg`, the YAML says `kind: csv, path:
foo.svg` and the CLI loader chokes.

Fix sketch in `notes/pending_work.md` item 3.

### Outline-mode title fallback

When uploaded files are in outline mode, the preview title says
`outline ↔ outline` (or `outline ↔ cap`) — uses literal string
instead of the uploaded filename. Cosmetic.

### Slider rebuild loses values on unit change

When you change the **Length units** dropdown, the slider cells
re-execute and recreate the inputs with new default values. Anything
you typed gets reset to the default for the new unit. Not technically
a bug, but surprising UX.

To preserve, you'd have to thread the previous values through a
state cell. Marimo can do this but the current code doesn't.

## CLI / mesher

### `dedup_consecutive` tolerance is unit-agnostic

`mesh/layout.py:dedup_consecutive(..., tol=1e-9)`. For µm-scale
geometry, tol = 1e-9 µm = 1 fm — almost no tolerance. Works for
all current cases because nothing has duplicate vertices that close,
but should scale with polyline size to be robust.

### Auto-scale guard threshold

`slice_2d.py` clamps `mesh_size` when it exceeds `4 × min_seg`. The
multiplier is somewhat arbitrary — picked to leave the existing
builtin examples unchanged (their `mesh_size / min_seg` ratio is
~2.4) while catching obviously-wrong inputs (>4×).

Should ideally be informed by gmsh's actual mesh-recovery
heuristics, not picked by hand.

### Hp = 0 corner case

For STEP-extracted right-half profiles where the body starts at
z=0 (not z=-Hp), the schema requires `Hp_um: gt=0` so the user has
to set `Hp_um: 0.001` (essentially zero). The mesher computes
`Hp = -half_z[0] = 0` from the polyline itself anyway, so
`Hp_um` field is ignored in this case.

Either drop the `gt=0` constraint OR have the mesher use `Hp_um`
directly when `half_z[0] != -Hp_um`. Current behavior is "works
but feels weird."

## Solver

### Newton convergence check is AND not OR

`dolfinx_contact.newton_solver.NewtonSolver` (third-party, lives in
the container, not in this repo) requires BOTH `r_abs < atol` AND
`r_rel < rtol`. Common convention is OR. This is why the pin/socket
sim stalls at step 6 even though `r_rel = 0.0906` is well below
`rtol = 0.2`.

Can't easily change without forking dolfinx_contact. Workaround:
set `atol` generously high.

### `gamma_scale` unit convention

The Nitsche penalty scales with `E_mpa * gamma_scale` — so the
absolute penalty depends on the material's stiffness. The defaults
were tuned for the µm-scale spheres with E = 200 GPa. For larger
parts you can usually leave γ alone (the E factor handles
non-dimensionalization) but for sharp barb geometry, smaller γ
gives better convergence.

Not a bug per se — just non-obvious. Documenting here so the next
maintainer doesn't go on a wild chase.

## STEP extractor

### Multi-face cross-sections may pick wrong face

`step_to_profile.py:extract_all_outlines` selects the plane
(z=0 or z=1) with the most face area, then walks all faces on that
plane. If your STEP has construction faces or sketch leftovers
mixed into the relevant plane, the extracted outline includes them.

Workarounds: clean the STEP (remove construction geometry) or
manually edit the resulting CSV/SVG.

### Single-face CSV from full outline

When `--full-outline` runs on a multi-lobe geometry (e.g. socket),
the SVG correctly emits one `<path>` per lobe. But the CSV only
emits the largest face's polygon — clients that need the full
geometry must use the SVG.

Symmetry of behavior would be nicer (CSV could emit multi-line-pair
sections separated by a blank line) but it's a documented limitation.

## Tests

### `test_marimo_inline_sync.py` is conservative

The inline-sync test only checks the SHAPES dict and the layout
helpers. The newer inline loaders (`load_csv_bytes`,
`load_csv_outline_bytes`, etc.) and the `outline_body_layout`
helper aren't sync-tested against package versions — because they
don't have package counterparts. If we add them to the package
later (for non-GUI consumers), add corresponding sync tests.

### No simulation correctness tests

The test suite covers polylines and configs but doesn't actually
run a solve. The first time something breaks the solver pipeline,
you'll only find out by running `./scripts/run_config.sh`. A
cheap smoke test would: run sphere_push.yaml's mesher (no solver),
check the mesh.msh has expected number of cells. Faster than
nothing.
