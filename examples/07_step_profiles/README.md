# Example 07 — STEP-derived profiles

CSV + SVG profile files produced from STEP slices via
`scripts/step_to_profile.py`. Once converted, these can be:

- Fed back through the CLI: `--profile-kind csv --profile-path examples/07_step_profiles/pin_slice.csv`
- Uploaded in the Marimo GUI (Bottom / Top profile source = `csv` or `svg`)

## Source files

- `pin_slice.csv`, `pin_slice.svg`  — extracted from `FEA_pin_slice.step`
  (1.0 mm slab, axis along STEP-x, axis at y=0)

## Re-generating

```bash
python scripts/step_to_profile.py /path/to/FEA_pin_slice.step \
  --csv examples/07_step_profiles/pin_slice.csv \
  --svg examples/07_step_profiles/pin_slice.svg \
  --units mm --extrusion-axis z --rev-axis x --axis-offset 0.0
```

The utility:

1. Loads the STEP via gmsh's OpenCASCADE kernel.
2. Picks the largest face flat in the extrusion direction (the 2D sketch silhouette).
3. Meshes its boundary as 1D segments and walks them topologically to
   recover ordered closed cycles (handles bodies with internal notches
   that fragment the boundary into multiple disjoint loops).
4. Picks the cycle with the largest |signed area| (the outer
   silhouette).
5. Splits at the symmetry-axis crossings and emits a right-half
   polyline in our (x_um, z_um) convention.

## Open: socket interpretation

`FEA_socket_slice.step` doesn't fit the single-closed-loop-around-an-axis
model — it has two separate lobes with a cavity between them. The
current pipeline expects a polyline that starts and ends on the
symmetry axis. Options:

- Extract one lobe only (treat its inner barbed surface as the contact
  profile, with the lobe's outer side and pillar shaft replaced by the
  standard pillar geometry).
- Extend the pipeline to accept full 2D outlines for non-axisymmetric
  bodies.

Need to decide which interpretation matches the physical intent before
generating a socket profile here.
