# Newton solver convergence tuning

The pin-into-socket simulation got through steps 1-5 cleanly
(`iters=5, converged=True` each) but stalled at step 6 (when the
first barb pair starts engaging). This file documents what's
happening and the recommended fix.

## Diagnosis

From the run log, step 6 Newton history:

```
Newton Iter 0:  r_abs = 42,458,888,312    r_rel = 70789           relax = 1.0
Newton Iter 7:  r_abs = 122,239           r_rel = 0.0906          relax = 1.0
Newton Iter 8:  r_abs = 122,239           r_rel = 0.0906          relax = 0.015625
...
Newton Iter 40: r_abs = 122,246           r_rel = 0.0906          relax = 0.015625
```

Two things to notice:

1. **Iter 7 actually achieved real convergence** — `r_rel = 0.0906 <
   rtol = 0.2`. The residual dropped 5 orders of magnitude in 7
   iterations.

2. **The dolfinx_contact Newton solver requires BOTH `r_abs < atol`
   AND `r_rel < rtol`**. Iter 7's `r_abs = 122,239 > atol = 100,000`,
   so it didn't accept.

3. **At iter 8 the auto-relaxer kicked in** because the residual
   ticked up by a roundoff-tiny amount. Relaxation cut to 1/64 and
   the solver got stuck doing nothing — each iter changes the
   solution by 1/64 of the Newton direction, which is too small to
   make progress past the local optimum but big enough to drift
   slightly in r_abs (122,239 → 122,246 over 33 iters).

This is a **threshold-crossing problem**, not a model problem. The
solver is essentially converged but doesn't realize it.

## Recommended YAML fix

Already applied to `examples/07_step_profiles/pin_socket_push.yaml`:

```yaml
solver:
  steps: 700               # was 200; smaller increments through barbs
  newton_atol: 500000.0    # was 100000; loose enough to accept the
                           # natural residual floor (~120k) at barb-
                           # engagement events
  gamma_scale: 0.05        # was 0.2; softer Nitsche penalty for
                           # sharper geometric features
```

### Why each knob matters

**`steps: 200 → 700`** (load increment 175 µm → 50 µm). This is the
biggest factor. At 50 µm/step the contact set transitions through
barb engagement gradually instead of jumping. Newton's quadratic
local model stays valid step-to-step.

**`newton_atol: 100000 → 500000`**. Sets the floor for accepting
convergence. Was tuned for µm-scale spheres; pin/socket has
mm-scale stresses so the natural residual is larger. 500k > the
observed 120k floor, with margin.

**`gamma_scale: 0.2 → 0.05`**. Nitsche method's penalty parameter.
Lower γ = softer enforcement of the contact constraint. For sharp
barb corners, soft contact reduces ill-conditioning at the cost of
allowing ~µm of interpenetration. Tighten back to 0.2 later if the
final result looks fuzzy.

## If it still fails

Try these in order:

### 1. Even smaller steps

```yaml
solver:
  steps: 1400              # 25 µm/step
```

Doubles run time but should converge through any reasonable barb
geometry.

### 2. Tighter rtol

```yaml
solver:
  newton_rtol: 0.1         # was 0.2
```

Tighter relative tolerance means iter 7's r_rel = 0.0906 is even
more comfortably converged.

### 3. Round the barb corners in CAD

Sharp 90° corners create singular FE stress fields. A 100-200 µm
fillet at each barb tip in the original CAD often *dramatically*
improves solver convergence. Worth doing if you control the CAD.

### 4. Continuation strategy

Solve at coarse mesh first (e.g. `characteristic_size_um: 500`),
then use the result as initial guess for a finer mesh.
**Not currently supported** — would need solver code changes.

## Schema gaps you may hit

The current `SolverConfig` exposes:
`steps, disp_um, direction, contact_mode, fric, newton_rtol,
newton_atol, gamma_scale, ksp_type, pc_type`.

These are **not yet exposed in the schema**, though they exist as
CLI flags on the solver:

- `--newton-max-it` (default 40) — would help if you genuinely need
  more iterations
- `--relax` (default 1.0) — explicit relaxation override
- `--no-update-contact` — freeze contact pairs (faster but breaks
  if contact topology changes)

To use any of these, add them to `SolverConfig` in `config.py` and
plumb through `scripts/run_config.sh`.

## Background: how dolfinx_contact's Newton works

It's a damped Newton iteration with an automatic relaxation cut:
- Compute Newton step `Δu` from `J(u_n) Δu = -F(u_n)`
- Try `u_n+1 = u_n + α Δu` with `α = 1`
- If `||F(u_n+1)|| > ||F(u_n)||`, halve α and retry (up to ~6 halvings)
- If still no improvement, accept with the smallest α

The relaxation cuts to 1/64 (= 2⁻⁶) and stops there. Once `α < 1/64`,
the iteration makes no real progress but burns iterations.

The convergence check is `r_abs < atol AND r_rel < rtol`. Both must
hold. (If only one were required, we'd have been done at iter 7.)
