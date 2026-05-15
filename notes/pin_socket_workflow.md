# Pin / socket simulation workflow

Step-by-step recipe for going from your STEP files to a simulation.
The pin and socket files used for development were:

- `/Users/mnakamura/Downloads/FEA_pin_slice.step` (12 mm wide × 39 mm tall)
- `/Users/mnakamura/Downloads/FEA_socket_slice.step` (two lobes, midline at y=1 mm)

Both are **extruded 2D sketches** (thin slabs in z, axis-of-revolution
along STEP-x). They turn out to be revolves, so the right-half
extraction is lossless.

## 1. Convert STEP → right-half CSV

```bash
# Pin: axisymmetric around STEP-y=0
python scripts/step_to_profile.py /path/to/FEA_pin_slice.step \
    --csv examples/07_step_profiles/pin_slice_rh.csv \
    --units mm --extrusion-axis z --rev-axis x --axis-offset 0.0 \
    --right-half --out-units um --mesh-size 0.5

# Socket: pick one lobe, set the axis to its outer edge (y=-6.5 for
# the lower lobe). The output is a "right-half" of the lobe, with
# the outer edge as the symmetry axis.
python scripts/step_to_profile.py /path/to/FEA_socket_slice.step \
    --csv examples/07_step_profiles/socket_slice_rh.csv \
    --units mm --extrusion-axis z --rev-axis x --axis-offset -6.5 \
    --right-half --out-units um --mesh-size 0.5
```

**Why `--mesh-size 0.5`** (mm)? Boundary segments end up ~100 µm
each — coarse enough that gmsh in the 2D mesher can sit ~50-100 µm
bulk triangles alongside them without "Impossible to recover edge"
warnings. Default `--mesh-size 0.1` (100 µm) gives ~50 µm boundary
segments which require <25 µm bulk triangles → millions of cells.

**Why `--right-half`?** The default is `--full-outline` for general
geometries. Right-half is required for the CLI mesher (which only
supports the axisymmetric-mirror model). Both pin and socket are
genuine revolves so right-half is lossless for them.

**Why `--out-units um`?** Keeps the CSV in µm, matching the
schema's `*_um` field convention. If you used mm, you'd need to set
`units: mm` in the YAML's profile block.

## 2. Write the YAML

Already done in `examples/07_step_profiles/pin_socket_push.yaml`
and `pin_socket_retention.yaml`. Key geometry:

- `pitch_um: 14000` — cell width slightly larger than pin width 11.6 mm
- `Hp_um: 0.001` — body has no separate pillar (z starts at 0, so
  `Hp = -half_z[0] = 0`). Schema requires `gt=0`, hence 0.001.
- `plate_um: 1000` — 1 mm thick base/top plates
- `initial_gap_um: -100` for push (slight engagement to start),
  `-30000` for retention (already deeply engaged)
- `mesh.characteristic_size_um: 200` — 200 µm triangles
- Tuned solver values: `steps: 700`, `newton_atol: 500000`,
  `gamma_scale: 0.05` (see `convergence_tuning.md` for why)

## 3. Run

```bash
./scripts/run_config.sh examples/07_step_profiles/pin_socket_push.yaml
```

Outputs land in `runs_step_profiles/pin_socket_push/`:
- `mesh/` — the gmsh mesh files
- `mesh_preview.png` — matplotlib view of the cell + boundary tags
- `results/force_displacement.csv` — Fz(disp) curve (if it converges)
- `results/contact_results.xdmf` — full displacement field per timestep
- `contact_animation_pv.gif` — animated visualization
- `run.log` — full stdout/stderr including Newton convergence info

## 4. What to check first

- **`mesh_preview.png`** — confirms the geometry is what you expect.
  Pin should be on the bottom with barbs facing up; socket lobe
  mirrored across `x=P/2` to fill the cell, sitting above pin with
  initial engagement.

- **`run.log`** — the Newton iterations. Each step line looks like:
  `step N/700: disp=X µm, Fz=Y µN/µm-depth, iters=K, converged=True`.
  If `converged=False`, the run aborts early.

- **`force_displacement.csv`** — engineering quantity of interest. If
  the run completed all 700 steps, this gives the push-in (or
  retention) force curve.

## 5. Common issues

### Mesher fails with "Expected triangle (2D cells) and line (facets)"

gmsh couldn't recover boundary edges. Causes:
- `mesh_size_um` >> smallest boundary segment → bulk triangles can't
  fit alongside boundary
- Boundary polyline has self-intersections after mirror-assembly

Fix: pass `--gmsh-terminal 1` to see gmsh's actual warnings, then
either reduce `mesh.characteristic_size_um` or regenerate the CSV
with a different `--mesh-size` value in `step_to_profile.py`.

### Mesher fails with "polyline endpoints not on axis" (CSV loader)

Means the CSV's first or last point isn't at `x≈0`. Right-half
profiles must start and end on the symmetry axis.

If your CSV is a full outline (closed loop), it can't be used
through the CLI yet — see `pending_work.md` for the outline-mode
mesher work needed.

### Newton solver doesn't converge

See `convergence_tuning.md` — the most likely fixes are smaller
`steps` increments and looser `newton_atol`.

## Alternative: full-outline workflow (preview only)

To visualize the actual STEP cross-section (not the right-half
approximation) in the Marimo GUI:

```bash
python scripts/step_to_profile.py /path/to/part.step \
    --csv part.csv --svg part.svg --units mm --out-units mm
```

(No `--right-half` flag — default is full-outline.) Upload the SVG
or CSV in the GUI; the auto-detect logic in the resolver will fall
back to outline mode. The preview renders correctly; **but the YAML
download will not run through the CLI pipeline yet.** Need the
outline-mode mesher work.
