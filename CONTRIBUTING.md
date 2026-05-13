# Contributing

Thanks for your interest in contributing.

## Reporting bugs and requesting features

Open an issue at
<https://github.com/nanosystemslab/mis-contact-fea/issues>. Include:

- What you were trying to do
- What happened (and what you expected)
- A minimal reproducer if possible (smallest mesh / shape / DISP combo
  that triggers the issue)
- Your platform (macOS / Linux), Docker version, Python version

For simulation convergence problems, please also attach:

- The `runs_*/<shape>/run.log` tail
- The relevant `force_displacement.csv` if it exists
- Your `INITIAL_GAP` / `DISP` / `STEPS` settings

## Development setup

```bash
git clone https://github.com/nanosystemslab/mis-contact-fea
cd mis-contact-fea
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install pre-commit
pre-commit install
```

The simulator itself runs inside the `dolfinx-contact:local` Docker image
(`docker/build_asimov_image.sh`) so you don't need a working dolfinx /
petsc / mpi4py stack on the host.

## Code style

- Python: [black](https://github.com/psf/black) (line length 100) +
  [isort](https://pycqa.github.io/isort/) (black profile). Both run
  automatically via pre-commit.
- Shell: [shellcheck](https://www.shellcheck.net/) at `--severity=warning`.
- Files: LF line endings, UTF-8, final newline (enforced by `.editorconfig`).

The pre-commit suite is intentionally relaxed — no docstring enforcement,
no strict typing, no coverage gate. Tighten as the project matures.

## Branching / PRs

- Branch from `main` with a descriptive name (`feat/...`, `fix/...`,
  `docs/...`).
- Keep PRs focused. If a change touches both the solver and the
  analysis pipeline, split it.
- Rebase rather than merge before opening the PR.
- All commits must be signed off (`git commit -s`).

## Reproducibility checklist

If your change affects simulation results:

- [ ] Re-run the affected shape(s) via `./scripts/run_final_sweep.sh`
- [ ] Re-run retention via `./scripts/run_retention_sweep.sh`
- [ ] Re-run analysis via `python3 src/analyze_forces.py`
- [ ] Update the peak table in `README.md` if the numbers change
- [ ] Note the change in the PR description with before/after values

## Adding a new shape

1. Add a generator function in `src/profiles.py` that returns
   `(x, z)` polyline arrays.
2. Register it in the `SHAPES` dict at the bottom of that file.
3. Calibrate effective diameter to 100 µm if the new shape has corner
   fillets (see the cone/doublecone bisection comments).
4. Add per-shape `INITIAL_GAP` / `DISP` / `STEPS` rows to both
   `scripts/run_final_sweep.sh` and `scripts/run_retention_sweep.sh`.
5. Use `python3 src/preview_positions.py` to visually tune the gap
   and displacement before committing to a full sweep.
