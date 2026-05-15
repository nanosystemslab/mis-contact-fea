# Session log — 2026-05-14 → 2026-05-15

Chronological record of what we built and why, with commit SHAs so you
can `git show` any of them to see the exact diff.

## Day 1 (2026-05-14)

### `e8c7d58` — Asymmetric top/bottom profile support to 2D mesher

The original mesher assumed both bodies use the same profile. Added
`top_profile` optional field to `SimulationConfig`. Cleanest place
for the asymmetry is `top_body_vertices(..., bottom_z_apex=None)` in
`mesh/layout.py` — when set, the top body's apex anchors to the
bottom's apex + `initial_gap` regardless of profile apex-height
mismatch.

Tests in `tests/test_asymmetric_layout.py` cover:
1. `bottom_z_apex=None` reproduces legacy single-profile behavior
   byte-for-byte (back-compat).
2. With sphere bottom + prolate top, the top apex still lands
   `initial_gap` above the bottom apex.

### `acccba6` — Asymmetric top/bottom pickers to Marimo GUI

Added **Asymmetric pair** checkbox + **Top shape** dropdown in the
GUI. Dual-preview renders the pair with the new
`bottom_z_apex` anchoring. YAML preview emits `top_profile:` only
when asymmetric is on (back-compat for symmetric configs).

### `0515992` — CSV/SVG upload + ZIP bundle to Marimo GUI

Added file-upload widgets for CSV/SVG, plus inline `load_csv_bytes`
and `load_svg_bytes` loaders that work in WASM (no filesystem
required). Download cell auto-detects uploads and bundles
`<config>.yaml + profiles/<file>` into a ZIP when present, falls
back to plain YAML otherwise. `svgpathtools` added to the marimo
notebook's inline deps.

### `4536e2c` — Mirror top profile picker beside bottom

User feedback: the top picker was only visible when "Asymmetric"
was checked. Restructured so the picker rendering mirrors the
bottom layout when asymmetric is on — top/bottom side-by-side.

## Day 2 (2026-05-15)

### `4433c54` — Per-body x-offset sliders; write profiles in input units

Two unrelated changes bundled:

1. **x-offset sliders**: `bottom_x_offset` and `top_x_offset` number
   inputs in Geometry-advanced. Negative = shift left, positive =
   shift right (in user units). Plate stays at cell width [0, P]
   so you can see the body drift off-center visually.

2. **Output units honor input units**: `scripts/step_to_profile.py`
   gained `--out-units` flag. Default = same as `--units`. So a
   `--units mm` STEP gets written as a `x_mm,z_mm` CSV by default,
   avoiding the double-1000× scaling bug where the GUI loader
   re-multiplied µm values by 1000 when the user picked "mm" in
   the units dropdown.

### `2a10016` — Force outline mode for closed-loop SVG/CSV uploads

The right-half CSV/SVG loaders trivially "accept" any closed-loop
polyline because the start vertex and end vertex are the same point
(both satisfy the axis-touch check, even when both are off-axis).
That made socket SVGs sneak through as right-half profiles → got
mirrored across the seam → produced bizarre socket-as-pin rendering.

Fix in the GUI resolver: after a right-half load, check
`start == end as a 2D point` — if yes, it's a full outline and we
re-load via the outline-mode loader instead.

### `611e0a7` — Flip-in-z checkboxes for bottom and top bodies

`bottom_flip_z` and `top_flip_z` mo.ui.checkbox under each profile
source. When checked, the body's polyline gets `(x, z)` → `(x[::-1],
-z[::-1])` before positioning (negate z, reverse to keep CCW). Works
in both right-half mode (the layout helper still works because
`half_z[0] = -Hp` semantics are preserved) and outline mode.

### Uncommitted work in working tree

After `611e0a7`, the rest of the day's work is in the working tree
(not yet pushed). Most consequential additions:

**`src/mis_contact_fea/mesh/slice_2d.py`**: auto-scale guard in
`build_mesh()`. Computes the min segment length across the assembled
bottom and top boundary polygons; if `args.mesh_size_um` is more
than 4× larger than that, prints a warning and clamps to
`min_seg × 2`. Threshold deliberately permissive — gmsh CAN handle
some excess; we only intervene when the ratio is clearly broken.

**`src/mis_contact_fea/profiles/from_csv.py`**: `_compact_flat_runs()`
helper. STEP-extracted right-half profiles trace the polygon's full
boundary including the FLAT z=min edge (the body's base). The
existing `bottom_body_vertices()` / `top_body_vertices()` assume the
polyline starts at the pillar shoulder and immediately transitions
to the body side — flat-bottom runs duplicate the plate's inner
surface and fold the polygon onto itself near the top body's plate.

Compact runs of vertices at exactly `z_min` or `z_max` down to just
the two endpoints of each run. Tolerance is `1e-9` (essentially
exact equality) so curved approaches near an analytic shape's apex
aren't touched. The CSV roundtrip tests still pass.

**`examples/07_step_profiles/`**:
- `pin_slice_rh.csv` / `socket_slice_rh.csv` — right-half profiles
  from the user's `FEA_pin_slice.step` / `FEA_socket_slice.step`,
  generated with `--right-half --out-units um --mesh-size 0.5`.
  (Note: `mesh-size 0.5 mm` was deliberate — gives ~100 µm boundary
  segments, large enough that gmsh can sit 50 µm bulk triangles
  alongside them.)
- `pin_socket_push.yaml` — push-in simulation config. After
  parameter tuning (steps=700, atol=500000, gamma=0.05), should
  converge through barb engagement. Currently set with those tuned
  values.
- `pin_socket_retention.yaml` — pull-out config, similar geometry,
  `contact_mode: closest` (rather than `raytracing`) since separation
  events need closest-point detection.

## Things to remember about the working state

1. **GitHub Pages is live** at
   <https://nanosystemslab.github.io/mis-contact-fea/>.
   The deploy workflow fires on push to main when
   `align_marimo.py` or `.github/workflows/deploy-gui.yml` changes.

2. **The `claude` contributor on GitHub got removed** when the
   repo was deleted + recreated. Going forward, my commit messages
   omit the `Co-Authored-By` trailer unless explicitly requested.

3. **Docker container** `mattnakamura/dolfinx-contact:v0.9.0` is
   the solver image. Image is locally cached (~5 GB).

4. **The pipeline is host-mesh + container-solve**:
   - host: `python -m mis_contact_fea.mesh.slice_2d`
   - container: `docker run mattnakamura/dolfinx-contact:v0.9.0 python -m mis_contact_fea.solver.contact_2d`
   - host: `python -m mis_contact_fea.postproc.animate`
   `scripts/run_config.sh` orchestrates all three.

5. **dolfinx is NOT in the host Python env** — only inside the
   container. Don't try to `import dolfinx` on the Mac.

6. **The CLI mesher / solver only handles right-half profiles.**
   The Marimo GUI handles outline mode for *preview*, but the
   downloaded YAML + the CLI pipeline still need the
   `bottom_body_vertices_outline` work outlined in
   `pending_work.md`.
