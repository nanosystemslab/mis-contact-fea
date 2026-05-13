# Contact FEA Runbook (basic_fea)

This file summarizes the workflow we established for building the contact-enabled container, meshing, and running the FEniCSx contact simulation locally and on KOA HPC.

## 1) Local Docker image (dolfinx_contact)

We added a Docker image that installs `dolfinx_contact` plus runtime deps.

### Build image
```bash
./docker/build_asimov_image.sh dolfinx/dolfinx:v0.9.0 dolfinx/dolfinx-asimov:0.9.0
```

### Test image
```bash
docker run --rm -it dolfinx/dolfinx-asimov:0.9.0 \
  python3 -c "import dolfinx_contact; print(dolfinx_contact.__file__)"
```

## 2) Meshing workflow (local or HPC)

Script: `src/build_contact_mesh.py`

Key features:
- Imports `*_base.STEP` + `*_structure.STEP`
- Aligns base bottom to z=0, centers on x/y
- Aligns structure on x/y and shifts in z
- `--auto-gap` (default) computes min vertical gap on contact surfaces and adjusts to requested clearance
- Writes: `mesh.msh`, `mesh.xdmf`, `facet_tags.xdmf`, `mesh_tags.json`
- Prints `Contact gap after alignment: ...` and surface selection counts

Example:
```bash
python3 src/build_contact_mesh.py \
  --base-step models/N1004_base.STEP \
  --structure-step models/N1004_structure.STEP \
  --out-dir mesh_out \
  --clearance-mm 1 \
  --mesh-size-mm 2 \
  --flip-xy
```

Optional alignment tweaks:
- `--x-offset-mm`, `--y-offset-mm` to shift structure in-plane
- `--z-offset-mm` to fine-tune vertical shift
- `--no-auto-gap` to disable auto gap alignment

## 3) Contact solve (local)

Script: `src/run_contact_sim.py`

Reads mesh from `mesh.msh` (preferred) or XDMF.
Uses dolfinx_contact Nitsche formulation with friction (default Coulomb).

Example:
```bash
python3 src/run_contact_sim.py \
  --mesh-dir mesh_out \
  --out-dir results_contact \
  --steps 20 \
  --disp-mm 1.2
```

Useful flags:
- `--fric 0.0` (start frictionless)
- `--gamma-scale 1.0` (lower penalty)
- `--relax 0.5` (Newton damping)
- `--newton-max-it 80`
- `--ksp-rtol`, `--ksp-atol`, `--ksp-max-it`

## 4) KOA HPC container (Singularity)

We follow the Optimization_Framework container style.

Files:
- `containers/dolfinx_v0.9.0_contact.def`
- `containers/build_contact_container.sh`
- `hpc/build_contact_container.slurm`

Build on KOA:
```bash
cd ~/basic_fea
sbatch hpc/build_contact_container.slurm
```

Output SIF:
```
~/basic_fea/containers/dolfinx_v0.9.0_contact.sif
```

The def file assumes the base SIF:
```
/home/mtdsn/Optimization_Framework/containers/dolfinx_v0.9.0.sif
```
Update `From:` if your base container is elsewhere.

## 5) KOA HPC run (Slurm)

Slurm script: `hpc/run_contact_sim.slurm`

Defaults:
- partition: `shared`
- time: `02:00:00`
- RANKS defaults to `SLURM_NTASKS`, but if `RANKS=1` it runs without MPI.
- email notifications enabled.

Submit:
```bash
cd ~/basic_fea
RANKS=1 sbatch -p shared -t 04:00:00 hpc/run_contact_sim.slurm
```

Common overrides:
```bash
SURFACE_TOL_MM=5 \
MESH_SIZE_MM=3 \
STEPS=10 \
DISP_MM=1.2 \
FRIC=0.0 \
GAMMA_SCALE=1.0 \
RELAX=0.5 \
NEWTON_MAX_IT=80 \
RANKS=1 \
sbatch -p shared -t 04:00:00 hpc/run_contact_sim.slurm
```

Monitoring:
```bash
squeue -u mtdsn
tail -f logs/contact_<JOBID>.log
tail -f logs/contact_<JOBID>.err
```

## 6) Known issues + fixes

### Missing facet tags
If logs say:
```
Missing facet tags: structure_contact, structure_top
```
Increase `SURFACE_TOL_MM` (e.g., 5) or rely on the built-in fallback tolerance.

### Newton not converging
If logs show `Newton solver did not converge`:
Use frictionless + lower gamma + damping + more Newton iterations:
```
FRIC=0.0 GAMMA_SCALE=1.0 RELAX=0.5 NEWTON_MAX_IT=80
```

### Job timeouts
Run on `shared` with longer time (`-t 04:00:00`) and/or reduce mesh size.

## 7) Files added/modified

- `src/build_contact_mesh.py` (auto-gap, robust tagging, offsets)
- `src/run_contact_sim.py` (robust tag checks, pinning, solver controls)
- `docker/Dockerfile.asimov` (adds dolfinx_contact deps)
- `docker/build_asimov_image.sh`
- `containers/dolfinx_v0.9.0_contact.def`
- `containers/build_contact_container.sh`
- `hpc/build_contact_container.slurm`
- `hpc/run_contact_sim.slurm`
- `usage.md`
