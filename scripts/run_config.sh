#!/usr/bin/env bash
# Run the full 2D contact-FEA pipeline (mesh + solve + animate) from a
# single YAML config.
#
# Usage:
#   ./scripts/run_config.sh examples/01_seven_builtins/sphere_push.yaml
#   ./scripts/run_config.sh examples/01_seven_builtins/cone_retention.yaml
#
# The output dir is read from the config's `output.out_dir` field.
# Inside that dir we write:
#     mesh/            (mesh.msh, mesh_tags.json, ...)
#     mesh_preview.png
#     results/         (force_displacement.csv, contact_results.xdmf)
#     contact_animation_pv.gif
#     run.log
#
# Implementation note: the HOST reads the YAML config (needs pydantic +
# pyyaml). The DOCKER container only sees CLI flags, so the dolfinx
# image doesn't need to know about YAML or pydantic.
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

CONFIG="${1:?usage: $0 path/to/config.yaml}"
if [ ! -f "$CONFIG" ]; then
    echo "config not found: $CONFIG" >&2
    exit 1
fi

# Container coordinates — overridable via env.
IMAGE="${IMAGE:-mattnakamura/dolfinx-contact:v0.9.0}"
CONTAINER_PYTHONPATH="/work/src:/usr/local/dolfinx-real/lib/python3.12/dist-packages:/usr/local/lib"

# Expand config into a flat set of shell-safe variables.
eval "$(python3 - "$CONFIG" <<'PY'
import shlex, sys
from pathlib import Path
sys.path.insert(0, "src")
from mis_contact_fea.config import SimulationConfig

cfg = SimulationConfig.from_yaml(sys.argv[1])
mapping = {
    "OUT_DIR":      cfg.output.out_dir,
    "MESH_DIR":     cfg.mesh_dir,
    "RESULTS_DIR":  cfg.results_dir,
    "STEPS":        cfg.solver.steps,
    "DISP":         cfg.solver.disp_um,
    "DIRECTION":    cfg.solver.direction,
    "CONTACT_MODE": cfg.solver.contact_mode,
    "FRIC":         cfg.solver.fric,
    "NEWTON_RTOL":  cfg.solver.newton_rtol,
    "NEWTON_ATOL":  cfg.solver.newton_atol,
    "GAMMA_SCALE":  cfg.solver.gamma_scale,
    "KSP_TYPE":     cfg.solver.ksp_type,
    "PC_TYPE":      cfg.solver.pc_type,
}
for k, v in mapping.items():
    print(f"{k}={shlex.quote(str(v))}")
PY
)"

MESH_PREVIEW="$OUT_DIR/mesh_preview.png"
GIF_OUT="$OUT_DIR/contact_animation_pv.gif"
LOG_FILE="$OUT_DIR/run.log"

mkdir -p "$OUT_DIR"

echo "==== Running config: $CONFIG ===="
echo "  out_dir = $OUT_DIR"
echo "  image   = $IMAGE"
echo

{
    echo "==> 1/4  Building mesh"
    python3 -m mis_contact_fea.mesh.slice_2d --config "$CONFIG"

    echo
    echo "==> 2/4  Rendering mesh preview"
    python3 -m mis_contact_fea.postproc.preview_mesh \
        --mesh-dir "$MESH_DIR" \
        --out "$MESH_PREVIEW"

    echo
    echo "==> 3/4  Solving contact problem"
    docker run --rm -v "$PWD:/work" -w /work -e PYTHONPATH="$CONTAINER_PYTHONPATH" "$IMAGE" \
        python3 -m mis_contact_fea.solver.contact_2d \
            --mesh-dir "$MESH_DIR" \
            --out-dir "$RESULTS_DIR" \
            --steps "$STEPS" \
            --disp-um "$DISP" \
            --fric "$FRIC" \
            --gamma-scale "$GAMMA_SCALE" \
            --newton-rtol "$NEWTON_RTOL" \
            --newton-atol "$NEWTON_ATOL" \
            --ksp-type "$KSP_TYPE" \
            --pc-type "$PC_TYPE" \
            --contact-mode "$CONTACT_MODE" \
            --direction "$DIRECTION"

    echo
    echo "==> 4/4  Animating"
    python3 -m mis_contact_fea.postproc.animate \
        --results-dir "$RESULTS_DIR" \
        --out "$GIF_OUT" \
        --fps "${FPS:-6}" \
        --warp-scale "${WARP_SCALE:-1}" \
        --clip-percentile "${CLIP_PCT:-99}" \
        --max-frames "${MAX_FRAMES:-60}" \
        --save-frames
} 2>&1 | tee "$LOG_FILE"

echo
echo "Done."
echo "  mesh:      $MESH_DIR/"
echo "  preview:   $MESH_PREVIEW"
echo "  results:   $RESULTS_DIR/"
echo "  animation: $GIF_OUT"
echo "  log:       $LOG_FILE"
