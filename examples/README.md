# Examples

Each directory is a self-contained workflow. Numbering reflects rough
order of complexity, not dependencies — pick whichever matches your use
case.

| # | example | when to use it |
|---|---------|----------------|
| **01** | [Seven built-in shapes](01_seven_builtins/) | Reproduce the published results table. 14 YAML configs (7 push + 7 retention) for the built-in `sphere / oblate / prolate / cone / cap / capcone / doublecone` lobes. |
| **02** | [Custom CSV profile](02_custom_csv_profile/) | You have a digitized `(x, z)` outline (from WebPlotDigitizer, measurement data, or a Python generator). Two-column CSV, any length unit. |
| **03** | [Custom SVG profile](03_from_svg_profile/) | You've sketched a profile in Inkscape / Illustrator / Figma. Loader parses `<path>` and converts to a polyline. |
| **04** | [Axisymmetric section of a STEP file](04_from_step_profile/) | You have a CAD model. Loader slices the body at y=0 and extracts the right-half profile — gives you 2D answers from 3D geometry without running a full 3D contact problem. |
| **05** | [3D contact demo](05_3d_demo/) | **Planned / not yet implemented.** True tetrahedral 3D contact; placeholder for the roadmap. |

## Quick reference

Run a single config:
```bash
./scripts/run_config.sh examples/01_seven_builtins/sphere_push.yaml
```

Run a glob:
```bash
./scripts/run_configs.sh examples/01_seven_builtins/*_push.yaml
```

Rebuild comparison plots from existing run output (no resimulation):
```bash
./scripts/build_comparison_plots.sh
```

Launch the interactive 2D alignment GUI to author new configs visually:
```bash
./scripts/launch_gui.sh
```

## How a config maps to a run

```
geometry      → mesh dimensions (pitch, pillar height, plate thickness, initial gap)
profile       → which lobe shape (builtin / csv / svg / step) and what length unit it's in
mesh          → characteristic gmsh element size
solver        → Newton tolerances, contact mode, total displacement, direction (push/retention)
output.out_dir → where mesh, results, animation, and run log land
```

See [`src/mis_contact_fea/config.py`](../src/mis_contact_fea/config.py)
for the full schema with docstrings.
