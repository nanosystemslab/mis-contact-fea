#!/usr/bin/env bash
# Run the full 2D contact pipeline for every shape in SHAPES, each into
# its own mesh / results / preview / gif directory.
#
# Override the per-shape solver settings via env vars (same names as
# run_2d_pipeline.sh):
#   PITCH INITIAL_GAP HP PLATE MESH_SIZE STEPS DISP FRIC GAMMA_SCALE
#   NEWTON_RTOL NEWTON_ATOL FPS CLIP_PCT WARP_SCALE IMAGE
#
# Default knobs match the "looks good at P=198" recipe.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Default sweep — pass SHAPES=... to override
SHAPES="${SHAPES:-sphere oblate prolate cone cap capcone doublecone}"

# Per-shape DISP and INITIAL_GAP defaults. Each shape has a different
# lobe height (z_apex), so the displacement needed to reach full
# interlock differs. Override any of these via env vars when invoking
# (e.g. DISP_prolate=160 ./run_all_shapes.sh) or set a global default
# DISP / INITIAL_GAP that applies to all shapes lacking an explicit
# override.
shape_disp() {
    case "$1" in
        sphere)     echo "${DISP_sphere:-${DISP:-100}}" ;;
        oblate)     echo "${DISP_oblate:-${DISP:-80}}" ;;     # Rz=25, shorter lobe
        prolate)    echo "${DISP_prolate:-${DISP:-160}}" ;;   # Rz=75, taller lobe
        cone)       echo "${DISP_cone:-${DISP:-130}}" ;;
        cap)        echo "${DISP_cap:-${DISP:-110}}" ;;
        capcone)    echo "${DISP_capcone:-${DISP:-160}}" ;;   # umbo on top: 75 µm tall
        doublecone) echo "${DISP_doublecone:-${DISP:-130}}" ;;
        *)          echo "${DISP:-100}" ;;
    esac
}
shape_steps() {
    case "$1" in
        sphere)     echo "${STEPS_sphere:-${STEPS:-50}}" ;;
        oblate)     echo "${STEPS_oblate:-${STEPS:-50}}" ;;
        prolate)    echo "${STEPS_prolate:-${STEPS:-50}}" ;;
        cone)       echo "${STEPS_cone:-${STEPS:-50}}" ;;
        cap)        echo "${STEPS_cap:-${STEPS:-50}}" ;;
        capcone)    echo "${STEPS_capcone:-${STEPS:-50}}" ;;
        doublecone) echo "${STEPS_doublecone:-${STEPS:-50}}" ;;
        *)          echo "${STEPS:-50}" ;;
    esac
}
shape_initial_gap() {
    case "$1" in
        sphere)     echo "${INITIAL_GAP_sphere:-${INITIAL_GAP:--50}}" ;;
        oblate)     echo "${INITIAL_GAP_oblate:-${INITIAL_GAP:--25}}" ;;
        prolate)    echo "${INITIAL_GAP_prolate:-${INITIAL_GAP:--75}}" ;;
        cone)       echo "${INITIAL_GAP_cone:-${INITIAL_GAP:--40}}" ;;
        cap)        echo "${INITIAL_GAP_cap:-${INITIAL_GAP:--40}}" ;;
        capcone)    echo "${INITIAL_GAP_capcone:-${INITIAL_GAP:--60}}" ;;
        doublecone) echo "${INITIAL_GAP_doublecone:-${INITIAL_GAP:--40}}" ;;
        *)          echo "${INITIAL_GAP:--50}" ;;
    esac
}

export STEPS="${STEPS:-50}"
export PITCH="${PITCH:-198}"
export HP="${HP:-220}"
export PLATE="${PLATE:-30}"
export MESH_SIZE="${MESH_SIZE:-4}"
export FRIC="${FRIC:-0}"
export GAMMA_SCALE="${GAMMA_SCALE:-1}"
export NEWTON_RTOL="${NEWTON_RTOL:-1e-5}"
export NEWTON_ATOL="${NEWTON_ATOL:-1e-3}"
export FPS="${FPS:-6}"
export CLIP_PCT="${CLIP_PCT:-99}"
export WARP_SCALE="${WARP_SCALE:-1}"
export IMAGE="${IMAGE:-dolfinx-contact:local}"

OUT_ROOT="${OUT_ROOT:-runs_all_shapes}"
mkdir -p "$OUT_ROOT"

SUMMARY="$OUT_ROOT/summary.txt"
{
    echo "Sweep started: $(date)"
    echo "Settings: PITCH=$PITCH STEPS=$STEPS"
    echo "          HP=$HP PLATE=$PLATE MESH_SIZE=$MESH_SIZE"
    echo "          FRIC=$FRIC GAMMA_SCALE=$GAMMA_SCALE"
    echo "Per-shape DISP / INITIAL_GAP:"
    for s in $SHAPES; do
        echo "  $s: DISP=$(shape_disp "$s")  INITIAL_GAP=$(shape_initial_gap "$s")"
    done
    echo
} | tee "$SUMMARY"

FAILED=()
SUCCEEDED=()
for s in $SHAPES; do
    echo
    echo "########################################################################"
    echo "##  shape = $s"
    echo "########################################################################"

    SHAPE_DIR="$OUT_ROOT/$s"
    mkdir -p "$SHAPE_DIR"

    export SHAPE="$s"
    export DISP="$(shape_disp "$s")"
    export INITIAL_GAP="$(shape_initial_gap "$s")"
    export STEPS="$(shape_steps "$s")"
    export MESH_DIR="$SHAPE_DIR/mesh"
    export RESULTS_DIR="$SHAPE_DIR/results"
    export MESH_PREVIEW="$SHAPE_DIR/mesh_preview.png"
    export GIF_OUT="$SHAPE_DIR/contact_animation_pv.gif"

    echo "    DISP=$DISP  INITIAL_GAP=$INITIAL_GAP"

    LOG_FILE="$SHAPE_DIR/run.log"
    if ./scripts/run_2d_pipeline.sh 2>&1 | tee "$LOG_FILE"; then
        SUCCEEDED+=("$s")
        if [ -f "$RESULTS_DIR/force_displacement.csv" ]; then
            cp "$RESULTS_DIR/force_displacement.csv" "$SHAPE_DIR/force_displacement.csv"
        fi
    else
        FAILED+=("$s")
        echo "!! $s FAILED, continuing to next shape"
    fi
done

echo
echo "########################################################################"
echo "Sweep complete"
echo "########################################################################"
{
    echo
    echo "Sweep finished: $(date)"
    echo "Succeeded (${#SUCCEEDED[@]}): ${SUCCEEDED[*]:-none}"
    echo "Failed    (${#FAILED[@]}): ${FAILED[*]:-none}"
    echo
    echo "Per-shape artifacts in $OUT_ROOT/<shape>/:"
    for s in $SHAPES; do
        echo "  $OUT_ROOT/$s/"
        echo "    mesh_preview.png       (mesh BCs + tags)"
        echo "    contact_animation_pv.gif"
        echo "    contact_animation_pv_frames/"
        echo "    force_displacement.csv"
        echo "    run.log"
    done
} | tee -a "$SUMMARY"
