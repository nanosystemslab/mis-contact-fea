# Working notes — mis-contact-fea

Snapshot of session state on **2026-05-15** so work can pick up on
another machine (or in a fresh AI session) without losing context.

## Where things stand

The project is a 2D plane-strain contact-FEA pipeline for mechanical
interlocking nanostructures (sphere, cone, cap, etc.). Over this
session we extended it from a builtin-shapes-only workflow into one
that can also handle:

- **Asymmetric pairs** (different bottom and top profiles)
- **Uploaded CSV / SVG profiles** in the Marimo GUI (CSV/SVG upload
  + ZIP bundle download)
- **STEP file extraction** to CSV/SVG via `scripts/step_to_profile.py`
- **Outline-mode bodies** in the GUI preview (full closed 2D
  cross-sections, not just right-half axisymmetric profiles)
- **Per-body x-offset and flip-in-z controls** in the GUI
- **Unit-aware geometry sliders** (mm / µm / cm / m) that rebuild
  with appropriate ranges when the units dropdown changes
- **Auto-scale guard in the mesher** that clamps `mesh_size` when
  the boundary polyline is way finer than the requested triangle size

We **proved the full pipeline runs end-to-end** for a pin-into-socket
simulation: STEP → CSV (via `step_to_profile.py --right-half`) → YAML
config → host mesher (gmsh) → Docker solver (dolfinx_contact) →
post-processing. The first 5 load-steps converged cleanly; step 6
(when the first barb pair starts engaging) fails to converge with the
default Newton tolerances — solver-tuning territory, not a structural
bug. See `convergence_tuning.md` for fix recommendations.

## Files in this directory

| File | What's in it |
|------|--------------|
| [`session_log.md`](session_log.md) | Chronological log of what we did this session, with the commits and reasoning |
| [`architecture.md`](architecture.md) | How the pipeline is structured (mesh / solver / GUI / docker) and which file does what |
| [`pin_socket_workflow.md`](pin_socket_workflow.md) | Step-by-step recipe for the user's actual workflow (STEP → CSV → sim) |
| [`convergence_tuning.md`](convergence_tuning.md) | Diagnosis of the step-6 convergence failure and recommended Newton-solver fixes |
| [`pending_work.md`](pending_work.md) | What's not done yet and rough effort estimates |
| [`known_bugs.md`](known_bugs.md) | Bugs we found but didn't fix, and where they live |

## Quick reference

```bash
# Local marimo GUI (same as the hosted Pages site, faster iteration)
marimo run src/mis_contact_fea/gui/align_marimo.py

# Convert a STEP file to a right-half CSV (works in CLI pipeline)
python scripts/step_to_profile.py /path/to/part.step \
    --csv out.csv --units mm --extrusion-axis z --rev-axis x \
    --axis-offset 0.0 --right-half --out-units um

# Convert a STEP file to a full-outline SVG (for visual inspection)
python scripts/step_to_profile.py /path/to/part.step \
    --csv out.csv --svg out.svg --units mm --out-units mm

# Run a full simulation (host mesher + Docker solver + animation)
./scripts/run_config.sh examples/01_seven_builtins/sphere_push.yaml

# Run the pin-into-socket example (needs convergence tuning — see notes)
./scripts/run_config.sh examples/07_step_profiles/pin_socket_push.yaml
```

## Uncommitted work in the working tree (as of snapshot)

These are local changes not yet pushed to GitHub:

- `src/mis_contact_fea/mesh/slice_2d.py` — auto-scale guard
- `src/mis_contact_fea/profiles/from_csv.py` — `_compact_flat_runs()`
- `examples/07_step_profiles/pin_slice_rh.csv` (new)
- `examples/07_step_profiles/socket_slice_rh.csv` (new)
- `examples/07_step_profiles/pin_socket_push.yaml` (new)
- `examples/07_step_profiles/pin_socket_retention.yaml` (new)
- `notes/` (this directory)

`runs_step_profiles/` is local run output and should stay gitignored.

Suggested commit on resume:
```bash
git add src/mis_contact_fea/mesh/slice_2d.py \
        src/mis_contact_fea/profiles/from_csv.py \
        examples/07_step_profiles/ notes/
git commit -m "Add pin/socket example + mesher auto-scale + CSV flat-run compaction"
```
