#!/usr/bin/env bash
# Iterative tweak + comparison-grid loop.
#
# Edits the per-shape overrides below, re-runs just those shapes, then
# rebuilds comparison_grid.png and opens it.
#
# Usage:
#   ./tweak_and_check.sh
#
# To re-run a different subset, override SHAPES on the command line:
#   SHAPES="oblate prolate cone" ./tweak_and_check.sh
#
# Per-shape overrides (edit these and re-run):
set -euo pipefail

PROJECT_DIR="/Users/mnakamura/Library/CloudStorage/OneDrive-TheAppliedResearchLabatUH/Documents/Nanosystems_Lab/basic_fea"
cd "$PROJECT_DIR"

# Which shapes to (re)run this iteration
SHAPES="${SHAPES:-oblate prolate}"

# Global tolerances + solver settings
export NEWTON_ATOL="${NEWTON_ATOL:-1.0}"
export NEWTON_RTOL="${NEWTON_RTOL:-5e-2}"
export GAMMA_SCALE="${GAMMA_SCALE:-0.3}"
export STEPS="${STEPS:-80}"

# Snapshot user-set overrides so we can re-apply them AFTER eval'ing the
# computed defaults (eval would otherwise clobber them).
USER_INITIAL_GAP_oblate="${INITIAL_GAP_oblate:-}"
USER_INITIAL_GAP_sphere="${INITIAL_GAP_sphere:-}"
USER_INITIAL_GAP_prolate="${INITIAL_GAP_prolate:-}"
USER_INITIAL_GAP_cone="${INITIAL_GAP_cone:-}"
USER_INITIAL_GAP_cap="${INITIAL_GAP_cap:-}"
USER_INITIAL_GAP_capcone="${INITIAL_GAP_capcone:-}"
USER_INITIAL_GAP_doublecone="${INITIAL_GAP_doublecone:-}"
USER_DISP_oblate="${DISP_oblate:-}"
USER_DISP_sphere="${DISP_sphere:-}"
USER_DISP_prolate="${DISP_prolate:-}"
USER_DISP_cone="${DISP_cone:-}"
USER_DISP_cap="${DISP_cap:-}"
USER_DISP_capcone="${DISP_capcone:-}"
USER_DISP_doublecone="${DISP_doublecone:-}"

echo "==== Step 1/3:  compute baseline INITIAL_GAP / DISP from sphere reference ===="
# Generate the equal-Euclidean-separation defaults (this exports
# INITIAL_GAP_<shape> and DISP_<shape> for all shapes).
eval $(python3 src/compute_initial_gaps.py \
        --pitch 198 \
        --baseline-shape sphere \
        --baseline-gap -50 \
        --export env)

# Now re-apply user overrides on top of the computed defaults. Built-in
# script defaults also live here — bump these between runs while tuning.
[ -n "$USER_INITIAL_GAP_oblate" ]    && export INITIAL_GAP_oblate="$USER_INITIAL_GAP_oblate"     || export INITIAL_GAP_oblate=10
[ -n "$USER_INITIAL_GAP_sphere" ]    && export INITIAL_GAP_sphere="$USER_INITIAL_GAP_sphere"
[ -n "$USER_INITIAL_GAP_prolate" ]   && export INITIAL_GAP_prolate="$USER_INITIAL_GAP_prolate"
[ -n "$USER_INITIAL_GAP_cone" ]      && export INITIAL_GAP_cone="$USER_INITIAL_GAP_cone"
[ -n "$USER_INITIAL_GAP_cap" ]       && export INITIAL_GAP_cap="$USER_INITIAL_GAP_cap"
[ -n "$USER_INITIAL_GAP_capcone" ]   && export INITIAL_GAP_capcone="$USER_INITIAL_GAP_capcone"
[ -n "$USER_INITIAL_GAP_doublecone" ] && export INITIAL_GAP_doublecone="$USER_INITIAL_GAP_doublecone"
[ -n "$USER_DISP_capcone" ]          && export DISP_capcone="$USER_DISP_capcone"               || export DISP_capcone=300
[ -n "$USER_DISP_prolate" ]          && export DISP_prolate="$USER_DISP_prolate"               || export DISP_prolate=340
[ -n "$USER_DISP_oblate" ]           && export DISP_oblate="$USER_DISP_oblate"
[ -n "$USER_DISP_sphere" ]           && export DISP_sphere="$USER_DISP_sphere"
[ -n "$USER_DISP_cone" ]             && export DISP_cone="$USER_DISP_cone"
[ -n "$USER_DISP_cap" ]              && export DISP_cap="$USER_DISP_cap"
[ -n "$USER_DISP_doublecone" ]       && export DISP_doublecone="$USER_DISP_doublecone"

echo "Per-shape overrides applied:"
for s in sphere oblate prolate cone cap capcone doublecone; do
    ig_var="INITIAL_GAP_$s"; d_var="DISP_$s"
    echo "  $s: INITIAL_GAP=${!ig_var}  DISP=${!d_var}"
done

echo
echo "==== Step 2/3:  run sweep for shapes: $SHAPES ===="
SHAPES="$SHAPES" ./run_all_shapes.sh

echo
echo "==== Step 3/3:  rebuild comparison_grid.png ===="
python3 src/compare_shapes.py \
    --root runs_all_shapes \
    --out comparison_grid.png

open comparison_grid.png
echo
echo "Done. comparison_grid.png opened."
