# Example 06 — asymmetric top/bottom profiles

Demonstrates the `top_profile:` field on `SimulationConfig`. The
bottom body uses one profile, the top uses another. Any of the
profile kinds (`builtin`, `csv`, `svg`, `step`) can appear on either
side independently.

`sphere_vs_cap.yaml` builds a unit cell with a sphere on the bottom
and a cap on the top. The top is positioned so its apex sits
`initial_gap_um` above the bottom apex, just like in the symmetric
case.

Build & run:

    python -m mis_contact_fea.mesh.slice_2d --config examples/06_asymmetric_pair/sphere_vs_cap.yaml
    python -m mis_contact_fea.solver.contact_2d --config examples/06_asymmetric_pair/sphere_vs_cap.yaml

Backward compatibility: configs without `top_profile` continue to
behave exactly as before (symmetric pair).
