"""3D contact-FEA solver — PLANNED, not yet implemented.

When implemented, this will be the 3D analogue of `solver.contact_2d`:
displacement-controlled Nitsche frictionless contact between two
tetrahedral bodies, running inside the `mattnakamura/dolfinx-contact`
image.

Expected scope at implementation time:

- Linear elasticity in 3D (E, ν as today).
- `dolfinx_contact.cpp.ContactMode.{ClosestPoint, Raytracing}`.
- Reaction force assembly on the top plate's surface (not edge as in 2D).
- Per-step CSV + XDMF time-series output, drop-in compatible with the
  existing `postproc.compare_forces` / `postproc.analyze` post-processing
  so the smoothing + R/I-ratio analysis works without changes.

Estimated effort: 2-4 weeks plus per-shape convergence tuning.
Compute: 10-100× the 2D plane-strain version; HPC strongly recommended.

See the "Roadmap" section in the top-level README.
"""
from __future__ import annotations


def main(*args, **kwargs):
    """Placeholder entry point."""
    raise NotImplementedError(
        "3D contact solver is not yet implemented. See "
        "https://github.com/nanosystemslab/mis-contact-fea/issues for the "
        "roadmap; mis_contact_fea.solver.contact_2d is the supported path."
    )


if __name__ == "__main__":
    raise SystemExit(
        "mis_contact_fea.solver.contact_3d is a roadmap placeholder. "
        "Use mis_contact_fea.solver.contact_2d for 2D plane-strain contact."
    )
