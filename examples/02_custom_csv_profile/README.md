# 02 — Custom CSV-defined profile

Demonstrates how to feed a custom right-half axisymmetric profile from
a CSV file, including non-µm units.

## What's here

```
teardrop_mm.csv   — full right-half polyline (pillar + lobe) in millimetres
config.yaml       — references the CSV with `profile.kind=csv, units=mm`
```

For the demo we use the same vertex set the built-in `sphere` generator
produces, but written in **millimetres** to show the unit conversion in
action. The pipeline result will be identical to the builtin sphere run.

## CSV format

A header row is auto-detected. Lines starting with `#` are treated as
comments. The polyline must:

- Start on the symmetry axis: first row has `x ≈ 0`
- End on the symmetry axis: last row has `x ≈ 0`
- Stay in the right half: all `x ≥ 0`

If your CSV contains **only the lobe outline** (and you want the
canonical pillar stitched to the bottom), add `prepend_pillar: true`
in the config's `profile:` block.

## Run it

```bash
./scripts/run_config.sh examples/02_custom_csv_profile/config.yaml
```

Output goes to `runs_custom/teardrop/` (as specified in the YAML).

## Adapting to your own shape

1. Digitize your lobe outline into two columns (x, z) — e.g. with
   [WebPlotDigitizer](https://apps.automeris.io/wpd/) from a literature
   figure, or read your own data.
2. Save as CSV. Pick whatever length unit is convenient (`um`, `mm`,
   `cm`, `m`) — set `profile.units` in the config to match.
3. Make sure the polyline starts and ends on the symmetry axis (x=0).
4. Either include the pillar in the CSV, or set
   `profile.prepend_pillar: true` and let the loader handle it.
