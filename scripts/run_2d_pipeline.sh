#!/usr/bin/env bash
# End-to-end 2D plane-strain interlock pipeline:
#   1) Build the unit-cell mesh
#   2) Render a quick mesh preview
#   3) Run the contact simulation in the dolfinx-contact container
#   4) Animate the result with pyvista
#
# Override any default by exporting the env var or passing on the CLI, e.g.
#   PITCH=160 DISP=140 ./run_2d_pipeline.sh
#
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Geometry
SHAPE="${SHAPE:-sphere}"
PITCH="${PITCH:-190}"
HP="${HP:-220}"
PLATE="${PLATE:-30}"
INITIAL_GAP="${INITIAL_GAP:--50}"
MESH_SIZE="${MESH_SIZE:-4}"

# Simulation
STEPS="${STEPS:-50}"
DISP="${DISP:-120}"
FRIC="${FRIC:-0}"
GAMMA_SCALE="${GAMMA_SCALE:-1}"
NEWTON_RTOL="${NEWTON_RTOL:-1e-5}"
NEWTON_ATOL="${NEWTON_ATOL:-1e-3}"
CONTACT_MODE="${CONTACT_MODE:-raytracing}"
DIRECTION="${DIRECTION:-down}"  # down = push-in, up = retention pull-out

# Output / docker
MESH_DIR="${MESH_DIR:-mesh_out_2d}"
RESULTS_DIR="${RESULTS_DIR:-results_2d_test}"
MESH_PREVIEW="${MESH_PREVIEW:-mesh_preview.png}"
GIF_OUT="${GIF_OUT:-contact_animation_pv.gif}"
IMAGE="${IMAGE:-mattnakamura/dolfinx-contact:v0.9.0}"

# Animation
FPS="${FPS:-6}"
CLIP_PCT="${CLIP_PCT:-99}"
WARP_SCALE="${WARP_SCALE:-1}"
MAX_FRAMES="${MAX_FRAMES:-60}"

cd "$PROJECT_DIR"

echo "==> 1/4  Building mesh  shape=$SHAPE pitch=$PITCH initial_gap=$INITIAL_GAP"
python3 -m mis_contact_fea.mesh.slice_2d \
    --shape "$SHAPE" \
    --pitch-um "$PITCH" \
    --Hp-um "$HP" \
    --plate-um "$PLATE" \
    --initial-gap-um "$INITIAL_GAP" \
    --mesh-size-um "$MESH_SIZE" \
    --out-dir "$MESH_DIR"

echo
echo "==> 2/4  Rendering mesh preview"
python3 -m mis_contact_fea.postproc.preview_mesh \
    --mesh-dir "$MESH_DIR" \
    --out "$MESH_PREVIEW"

echo
echo "==> 3/4  Solving contact problem  steps=$STEPS disp=$DISP gamma=$GAMMA_SCALE"
# The container ships PYTHONPATH=/usr/local/dolfinx-real/...:/usr/local/lib
# baked into the image (that's how dolfinx itself is on the path). We
# need to prepend /work/src so our package is also importable. Setting
# `-e PYTHONPATH=...` with `docker run` REPLACES the env var (no shell
# expansion happens), so we must include the full chain explicitly.
CONTAINER_PYTHONPATH="/work/src:/usr/local/dolfinx-real/lib/python3.12/dist-packages:/usr/local/lib"
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
        --ksp-type preonly \
        --pc-type lu \
        --contact-mode "$CONTACT_MODE" \
        --direction "$DIRECTION"

echo
echo "==> 4/4  Animating  fps=$FPS"
python3 -m mis_contact_fea.postproc.animate \
    --results-dir "$RESULTS_DIR" \
    --out "$GIF_OUT" \
    --fps "$FPS" \
    --warp-scale "$WARP_SCALE" \
    --clip-percentile "$CLIP_PCT" \
    --max-frames "$MAX_FRAMES" \
    --save-frames

echo
echo "Done. Artifacts:"
echo "  $PROJECT_DIR/$MESH_DIR/"
echo "  $PROJECT_DIR/$RESULTS_DIR/"
echo "  $PROJECT_DIR/$MESH_PREVIEW"
echo "  $PROJECT_DIR/$GIF_OUT"
