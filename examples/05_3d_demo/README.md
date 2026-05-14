# 05 — 3D contact demo (planned)

**Status:** roadmap / not yet implemented.

The 2D plane-strain pipeline (examples 01–04) is the supported path
today. This directory is a placeholder for the future 3D-tetrahedral
contact pipeline.

## Why 3D matters

The 2D plane-strain assumption is exact for an infinitely-deep
periodic slab, but real interlocking patches have finite extent and
edge effects. For patches narrower than ~10× the lobe diameter, or
geometries with non-axisymmetric features (e.g. anisotropic etching),
a true 3D run is the right call.

## What's coming

The 3D pipeline will mirror the 2D one:

- `mis_contact_fea.mesh.solid_3d` — gmsh tetrahedralisation, fed from
  built-in revolved profiles or directly from a STEP body.
- `mis_contact_fea.solver.contact_3d` — `dolfinx_contact` 3D contact,
  running inside the same `mattnakamura/dolfinx-contact` container.
- `mis_contact_fea.gui.align_3d` — PyVista / trame interactive viewer
  for positioning bodies in 3D before committing to a long run.

Configuration will use the same `SimulationConfig` schema (likely with
a new `geometry.dim: 3` field).

## Why it's not done yet

The 2D pipeline already takes 2-10 minutes per shape; 3D will be
10-100× more expensive (memory and wall time). It needs HPC, and
the Newton tolerances we tuned for 2D will need to be re-tuned for
the larger contact pair count. Best estimate: 2-4 weeks of focused
work + ongoing per-shape tuning.

## Track / contribute

Open or follow an issue on
[github.com/nanosystemslab/mis-contact-fea](https://github.com/nanosystemslab/mis-contact-fea/issues).
PRs that implement the 3D pipeline are welcome — see
`src/mis_contact_fea/{mesh/solid_3d.py,solver/contact_3d.py,gui/align_3d.py}`
for the expected entry points.

## In the meantime: 2D from a 3D STEP file

If you have a CAD model and just need a 2D answer fast, use the
**STEP axisymmetric section** path from `examples/04_from_step_profile/`.
That gives you ~80% of the value of full 3D at ~5% of the cost.
