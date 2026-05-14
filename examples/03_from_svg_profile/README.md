# 03 — Custom SVG-defined profile

Use a vector-drawn outline from Inkscape / Illustrator / Figma to seed
a simulation. Useful when you have a hand-sketched or
literature-traced shape that's hard to express analytically.

## What you need

1. An SVG file containing **one `<path>` element** that traces the
   right-half outline of your lobe (and optionally the pillar). The
   first and last vertices of the path must sit on the symmetry axis
   (the y-axis in screen coordinates).
2. A length unit declaration — the SVG numbers are dimensionless, so
   you tell the loader what they represent via `profile.units`.

## Drawing tips (Inkscape)

- Set document units (File → Document Properties) to your target
  scale. If you set them to **millimetres**, use `profile.units: mm`.
- Draw your right-half profile. **Snap the start and end nodes to the
  y-axis** (x=0).
- Avoid sub-paths — the loader takes only the first one. Either
  combine shapes (Path → Combine) or simplify before exporting.
- Save as "Plain SVG" (File → Save As → Plain SVG).

## Example config

```yaml
geometry:
  pitch_um: 198.0
  initial_gap_um: -55.0
profile:
  kind: svg
  path: my_outline.svg
  units: mm
  prepend_pillar: true    # stitch the standard pillar if your SVG is lobe-only
mesh:
  characteristic_size_um: 4.0
solver:
  steps: 200
  disp_um: 50.0
  direction: down
  contact_mode: raytracing
output:
  out_dir: runs_custom/my_svg_shape
```

## Run it

```bash
pip install svgpathtools           # one-time dep
./scripts/run_config.sh examples/03_from_svg_profile/config.yaml
```

## Common issues

- **"polyline endpoints not on axis after parsing"** — your SVG path
  doesn't start and end on x=0 in the SVG's coordinate space. Open it
  in Inkscape, hit the leftmost vertex, and verify its X is 0 (or
  snap it there).
- **Strange scale in the result** — `profile.units` doesn't match what
  one Inkscape document unit represents. Check Document Properties →
  Display units.
- **Multiple subpaths** — the loader picks only `path_index=0`. Combine
  paths in Inkscape (Path → Combine, Ctrl+K) before exporting.
