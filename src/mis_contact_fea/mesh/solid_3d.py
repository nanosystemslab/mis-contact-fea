"""3D solid mesh generation — PLANNED, not yet implemented.

When implemented, this module will produce a 3D tetrahedral mesh of the
two interlocking bodies starting from one of:

- A built-in shape (revolved around the z-axis from a `profiles.SHAPES` entry).
- A STEP / STP body loaded via `gmsh.merge` + `getEntities`.
- An OBJ / STL surface mesh remeshed by gmsh.

The expected interface mirrors `mesh.slice_2d`:

    from mis_contact_fea.config import SimulationConfig
    from mis_contact_fea.mesh.solid_3d import build_mesh
    cfg = SimulationConfig.from_yaml("examples/05_3d_demo/config.yaml")
    build_mesh(cfg)

For details on the expected configuration surface, see the
"3D contact" section of the top-level README.
"""
from __future__ import annotations


def build_mesh(*args, **kwargs):
    """Placeholder. Raises until 3D meshing lands."""
    raise NotImplementedError(
        "3D solid meshing is not yet implemented. See "
        "https://github.com/nanosystemslab/mis-contact-fea/issues for the "
        "roadmap; the 2D pipeline (mesh.slice_2d) is the supported path "
        "for now."
    )


if __name__ == "__main__":
    raise SystemExit(
        "mis_contact_fea.mesh.solid_3d is a roadmap placeholder. "
        "Use mis_contact_fea.mesh.slice_2d for 2D plane-strain meshing."
    )
