#!/usr/bin/env bash
# Retention sweep: starts bodies in their FINAL push-in configuration
# (engaged state), then pulls the top plate UP to measure pull-out
# force. The displacement direction is reversed (+z instead of -z).
#
# For each shape, the retention initial-gap is computed automatically
# from the push parameters:
#     retention_INITIAL_GAP = push_INITIAL_GAP - push_DISP
# So if you ran a push with INITIAL_GAP=-55, DISP=85 (sphere), the
# retention starts at INITIAL_GAP=-140 — i.e., bodies start in the
# state they were left at the end of the push sim.
#
# Output goes to runs_retention/<shape>/ so it doesn't overwrite the
# push results in runs_all_shapes/.
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# ============================================================
# Push parameters per shape (must match what you pushed with).
# Edit these to match your push sweep.
# ============================================================
PUSH_INITIAL_GAP_sphere=-55
PUSH_INITIAL_GAP_oblate=-30
PUSH_INITIAL_GAP_prolate=-80
PUSH_INITIAL_GAP_cone=-75
PUSH_INITIAL_GAP_cap=-65
PUSH_INITIAL_GAP_capcone=-110
PUSH_INITIAL_GAP_doublecone=-80

# Retention DISP is shorter than push DISP — we only need to travel just
# past the disengagement point (peak unlock force + ~3 µm). Cap/capcone
# keep their original DISP because they hadn't disengaged yet in the
# previous run. STEPS chosen for ~0.2 µm/step uniform resolution.
PUSH_DISP_sphere=50
PUSH_DISP_oblate=25
PUSH_DISP_prolate=55
PUSH_DISP_cone=18
PUSH_DISP_cap=40
PUSH_DISP_capcone=45
PUSH_DISP_doublecone=28

PUSH_STEPS_sphere=250
PUSH_STEPS_oblate=125
PUSH_STEPS_prolate=275
PUSH_STEPS_cone=90
PUSH_STEPS_cap=200
PUSH_STEPS_capcone=225
PUSH_STEPS_doublecone=140

# Note: PUSH_INITIAL_GAP above is the PUSH starting gap. Retention
# initial gap is computed below as (push_IG - push_DISP_full_push). Use
# the full push DISP for the IG calculation, not the trimmed retention
# DISP, because the bodies start at the end-of-push position.
FULL_PUSH_DISP_sphere=85
FULL_PUSH_DISP_oblate=40
FULL_PUSH_DISP_prolate=120
FULL_PUSH_DISP_cone=30.9
FULL_PUSH_DISP_cap=40
FULL_PUSH_DISP_capcone=45
FULL_PUSH_DISP_doublecone=40
# ============================================================

# Solver tolerances — match the push sweep. The tight defaults in
# run_all_shapes.sh cause Newton to fail at the noise floor once the
# bodies disengage and there's no real residual to drive down.
export NEWTON_ATOL="${NEWTON_ATOL:-1000}"
export NEWTON_RTOL="${NEWTON_RTOL:-0.2}"
export GAMMA_SCALE="${GAMMA_SCALE:-0.2}"

export OUT_ROOT="${OUT_ROOT:-runs_retention}"
export DIRECTION="${DIRECTION:-up}"           # retention = pull-out
export CONTACT_MODE="${CONTACT_MODE:-closest}"
export DISP_SCALE="${DISP_SCALE:-1.0}"

SHAPES="${SHAPES:-sphere oblate prolate cone cap capcone doublecone}"

echo "==== Retention sweep settings ===="
echo "  OUT_ROOT      = $OUT_ROOT"
echo "  DIRECTION     = $DIRECTION   (top plate moves +z, pulling bodies apart)"
echo "  CONTACT_MODE  = $CONTACT_MODE"
echo "  SHAPES        = $SHAPES"
echo
echo "  Per-shape (retention IG = push_IG - full_push_DISP; retention DISP trimmed past disengagement):"
for s in $SHAPES; do
    push_ig_var="PUSH_INITIAL_GAP_$s";   push_ig="${!push_ig_var}"
    full_d_var="FULL_PUSH_DISP_$s";      full_d="${!full_d_var}"
    ret_d_var="PUSH_DISP_$s";            ret_d="${!ret_d_var}"
    ret_n_var="PUSH_STEPS_$s";           ret_n="${!ret_n_var}"
    retention_ig=$(python3 -c "print(round(${push_ig} - ${full_d}, 3))")
    export "INITIAL_GAP_$s"="$retention_ig"
    export "DISP_$s"="$ret_d"
    export "STEPS_$s"="$ret_n"
    echo "    $s  INITIAL_GAP=$retention_ig  DISP=$ret_d  STEPS=$ret_n"
done
echo

# Use the same per-shape machinery as the push sweep, just with
# OUT_ROOT and DIRECTION overridden via env.
export SHAPES
./scripts/run_all_shapes.sh

echo
n_shapes=$(echo "$SHAPES" | wc -w | tr -d ' ')
if [ "$n_shapes" -le 1 ]; then
    only_shape=$(echo "$SHAPES" | awk '{print $1}')
    echo "==== Single-shape retention run ($only_shape) — skipping comparison plots ===="
    echo "Per-shape artifacts:"
    echo "  $OUT_ROOT/$only_shape/mesh_preview.png"
    echo "  $OUT_ROOT/$only_shape/contact_animation_pv.gif"
    echo "  $OUT_ROOT/$only_shape/results/force_displacement.csv"
    open "$OUT_ROOT/$only_shape/mesh_preview.png" \
         "$OUT_ROOT/$only_shape/contact_animation_pv.gif" 2>/dev/null || true
else
    echo "==== Building retention force/comparison plots ===="
    python3 -m mis_contact_fea.postproc.compare_forces --abs --mode contact-aligned --contact-threshold 1 \
        --root "$OUT_ROOT" --out force_comparison_retention.png || true
    python3 -m mis_contact_fea.postproc.compare_shapes --root "$OUT_ROOT" \
        --out comparison_grid_retention.png || true

    open force_comparison_retention.png comparison_grid_retention.png 2>/dev/null || true
    echo "Done. Retention artifacts in $OUT_ROOT/, plots in force_comparison_retention.png and comparison_grid_retention.png"
fi
