#!/usr/bin/env bash
# Push sweep across all 7 shapes using the per-shape INITIAL_GAP / DISP /
# STEPS values dialed in via preview_positions.py. Step counts give
# ~0.2 µm/step uniform resolution.
#
# Override anything via env vars when invoking, e.g.:
#   SHAPES="cone doublecone" ./run_final_sweep.sh
#   DISP_sphere=90 STEPS_sphere=450 ./run_final_sweep.sh
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# ============================================================
# Per-shape push settings (matched to preview_positions.py).
# Edit these values to retune; run_retention_sweep.sh reads its own copy.
# ============================================================
PUSH_INITIAL_GAP_sphere=-55
PUSH_INITIAL_GAP_oblate=-30
PUSH_INITIAL_GAP_prolate=-80
PUSH_INITIAL_GAP_cone=-75
PUSH_INITIAL_GAP_cap=-65
PUSH_INITIAL_GAP_capcone=-110
PUSH_INITIAL_GAP_doublecone=-80

PUSH_DISP_sphere=85
PUSH_DISP_oblate=40
PUSH_DISP_prolate=120
PUSH_DISP_cone=30.9
PUSH_DISP_cap=40
PUSH_DISP_capcone=45
PUSH_DISP_doublecone=40

PUSH_STEPS_sphere=425
PUSH_STEPS_oblate=200
PUSH_STEPS_prolate=600
PUSH_STEPS_cone=155
PUSH_STEPS_cap=200
PUSH_STEPS_capcone=225
PUSH_STEPS_doublecone=200
# ============================================================

# Solver tolerances tuned for the harder contact cases (cone/doublecone).
export NEWTON_ATOL="${NEWTON_ATOL:-1000}"
export NEWTON_RTOL="${NEWTON_RTOL:-0.2}"
export GAMMA_SCALE="${GAMMA_SCALE:-0.2}"

export OUT_ROOT="${OUT_ROOT:-runs_all_shapes}"
export DIRECTION="${DIRECTION:-down}"           # push-in
export CONTACT_MODE="${CONTACT_MODE:-raytracing}"

SHAPES="${SHAPES:-sphere oblate prolate cone cap capcone doublecone}"

echo "==== Push sweep settings ===="
echo "  OUT_ROOT      = $OUT_ROOT"
echo "  DIRECTION     = $DIRECTION   (top plate moves -z, pushing in)"
echo "  CONTACT_MODE  = $CONTACT_MODE"
echo "  NEWTON_ATOL   = $NEWTON_ATOL"
echo "  NEWTON_RTOL   = $NEWTON_RTOL"
echo "  GAMMA_SCALE   = $GAMMA_SCALE"
echo "  SHAPES        = $SHAPES"
echo
echo "  Per-shape settings:"
for s in $SHAPES; do
    # Honor env-var override if the user passed one in.
    ig_user_var="INITIAL_GAP_$s";  ig_user="${!ig_user_var-}"
    d_user_var="DISP_$s";          d_user="${!d_user_var-}"
    n_user_var="STEPS_$s";         n_user="${!n_user_var-}"
    ig_def_var="PUSH_INITIAL_GAP_$s";  ig_def="${!ig_def_var}"
    d_def_var="PUSH_DISP_$s";          d_def="${!d_def_var}"
    n_def_var="PUSH_STEPS_$s";         n_def="${!n_def_var}"
    ig="${ig_user:-$ig_def}"
    d="${d_user:-$d_def}"
    n="${n_user:-$n_def}"
    export "INITIAL_GAP_$s"="$ig"
    export "DISP_$s"="$d"
    export "STEPS_$s"="$n"
    echo "    $s  INITIAL_GAP=$ig  DISP=$d  STEPS=$n"
done
echo

export SHAPES
./scripts/run_all_shapes.sh

echo
n_shapes=$(echo "$SHAPES" | wc -w | tr -d ' ')
if [ "$n_shapes" -le 1 ]; then
    # Single shape — open that shape's own artifacts; comparison plots
    # would just be one curve in a "grid" so we skip them.
    only_shape=$(echo "$SHAPES" | awk '{print $1}')
    echo "==== Single-shape run ($only_shape) — skipping comparison plots ===="
    echo "Per-shape artifacts:"
    echo "  $OUT_ROOT/$only_shape/mesh_preview.png"
    echo "  $OUT_ROOT/$only_shape/contact_animation_pv.gif"
    echo "  $OUT_ROOT/$only_shape/results/force_displacement.csv"
    open "$OUT_ROOT/$only_shape/mesh_preview.png" \
         "$OUT_ROOT/$only_shape/contact_animation_pv.gif" 2>/dev/null || true
else
    echo "==== Building force/comparison plots ===="
    python3 -m mis_contact_fea.postproc.compare_forces --abs --mode contact-aligned --contact-threshold 1 \
        --root "$OUT_ROOT" --out force_comparison_aligned.png || true
    python3 -m mis_contact_fea.postproc.compare_shapes --root "$OUT_ROOT" --out comparison_grid.png || true

    open force_comparison_aligned.png comparison_grid.png 2>/dev/null || true
    echo "Done. Push artifacts in $OUT_ROOT/, plots in force_comparison_aligned.png and comparison_grid.png"
fi
