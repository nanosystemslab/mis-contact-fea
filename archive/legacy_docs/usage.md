# Usage

## Run with the contact-enabled image

1) Mesh

```bash
docker run --rm -it -v "$PWD":/work -w /work dolfinx/dolfinx-asimov:0.9.0 \
  python3 src/build_contact_mesh.py \
  --base-step models/N1004_base.STEP \
  --structure-step models/N1004_structure.STEP \
  --out-dir mesh_out \
  --clearance-mm 30 \
  --mesh-size-mm 2
```

2) Solve

```bash
docker run --rm -it -v "$PWD":/work -w /work dolfinx/dolfinx-asimov:0.9.0 \
  python3 src/run_contact_sim.py \
  --mesh-dir mesh_out \
  --out-dir results_contact \
  --steps 10
```
