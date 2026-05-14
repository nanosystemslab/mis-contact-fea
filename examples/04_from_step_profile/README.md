# 04 — STEP file axisymmetric section

Cut a 2D right-half profile out of a 3D STEP body and feed it to the
plane-strain solver. Useful when you have a CAD model of the
structure you want to simulate.

This is the "shortcut" path mentioned in the project README — you get
arbitrary 3D geometry as input *without* paying the 10–100× cost of a
true 3D contact run. The catch: only works for axisymmetric bodies
(rotation-of-revolution geometry).

## Assumptions

- Body is axisymmetric around the **z-axis** (axis of revolution = z).
- Body is centered at x=0 in the STEP file.
- The STEP file contains a single body (the loader picks the first
  volume if there are several).

## Example config

```yaml
geometry:
  pitch_um: 198.0
  initial_gap_um: -55.0
profile:
  kind: step
  path: my_lobe.step
  units: mm           # whatever your STEP file is authored in
mesh:
  characteristic_size_um: 4.0
solver:
  steps: 200
  disp_um: 50.0
  direction: down
  contact_mode: raytracing
output:
  out_dir: runs_custom/my_step_shape
```

## Run it

```bash
./scripts/run_config.sh examples/04_from_step_profile/config.yaml
```

## How the section is built

The loader:

1. Imports the STEP body via gmsh's OpenCASCADE kernel.
2. Builds a thin slab at y=0 (the symmetry plane).
3. Intersects the body with the slab — this leaves a 2D face on y=0.
4. Walks the face's boundary curves and pulls vertex coordinates.
5. Deduplicates, sorts by z, and snaps the polyline endpoints to x=0.

The result is the same `(x, z)` polyline format the built-in
generators produce, so it slots into the rest of the pipeline
unchanged.

## Common issues

- **"STEP section produced no points"** — the body isn't actually
  centered on z, or the y-axis pass doesn't intersect it. Open in
  FreeCAD / Fusion to verify the orientation.
- **Multiple bodies** — the loader uses only the first. Combine bodies
  in CAD before exporting, or pre-cut to a single body.
- **Wrong units / scale** — STEP files do declare their unit but it's
  often wrong. Verify with `profile.units` matching what the body's
  numbers should mean (mm is the most common authoring unit in CAD).
