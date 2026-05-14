# 01 — Seven builtin shapes

14 YAML configs covering the original sweep: each of the seven built-in
lobe shapes (sphere, oblate, prolate, cone, cap, capcone, doublecone)
in both **insertion** (`_push.yaml`) and **retention** (`_retention.yaml`)
directions.

These are the exact parameters used to produce the results table in the
top-level `README.md`.

## Run a single config

```bash
./scripts/run_config.sh examples/01_seven_builtins/sphere_push.yaml
```

Writes mesh, results, animation, and run log to
`runs_all_shapes/sphere/` (or `runs_retention/sphere/` for retention
configs — the output directory is set inside each YAML's `output.out_dir`).

## Run a batch

```bash
# all 7 push configs
./scripts/run_configs.sh examples/01_seven_builtins/*_push.yaml

# all 7 retention configs (start each shape in its end-of-push state)
./scripts/run_configs.sh examples/01_seven_builtins/*_retention.yaml

# everything
./scripts/run_configs.sh examples/01_seven_builtins/*.yaml
```

## Build the comparative plots

After all 14 runs:

```bash
./scripts/build_comparison_plots.sh
```

Produces `force_smoothed_overlay.png`, `force_ratio.png`, `force_peaks.csv`
at the project root.

## Per-shape settings

| shape       | push IG | push DISP | push STEPS | retention IG | retention DISP | retention STEPS |
|-------------|---------|-----------|------------|--------------|----------------|-----------------|
| sphere      | −55     | 85        | 425        | −140         | 50             | 250             |
| oblate      | −30     | 40        | 200        | −70          | 25             | 125             |
| prolate     | −80     | 87        | 435        | −200 *       | 55             | 275             |
| cone        | −75     | 30.9      | 155        | −105.9       | 18             | 90              |
| cap         | −65     | 40        | 200        | −105         | 40             | 200             |
| capcone     | −110    | 45        | 225        | −155         | 45             | 225             |
| doublecone  | −80     | 40        | 200        | −120         | 28             | 140             |

\* Prolate retention starts at IG=−200 (which corresponds to a "full"
push of 120 µm from the IG=−80 starting position), even though the
actual push run only reached 87 µm before Newton chatter set in. This
mirrors the previous shell-script behavior; you can lower the
`initial_gap_um` value in `prolate_retention.yaml` to −167 if you want
the retention to start exactly where the push ended.

## Schema reference

See `src/mis_contact_fea/config.py` for the pydantic schema and full
docstrings. The validator catches typos and out-of-range values, so an
edit like `direction: dwon` will fail immediately rather than produce
silent garbage.
