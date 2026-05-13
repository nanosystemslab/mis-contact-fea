#!/usr/bin/env bash
# Revert sweep: roll back the two recent changes (Raytracing contact
# mode + early-stop logic) that may have broken the previously-good
# results, while keeping the apex fillets and per-shape bounds.
#
# Use this to get back to a known-working state. After it runs, compare
# the regenerated GIFs and force plots against what you remember from
# the working run.
#
# Skips prolate and cap (still under investigation).
set -euo pipefail

PROJECT_DIR="/Users/mnakamura/Library/CloudStorage/OneDrive-TheAppliedResearchLabatUH/Documents/Nanosystems_Lab/basic_fea"
cd "$PROJECT_DIR"

# --- knob that reverts the recent contact-mode change ---
export CONTACT_MODE="${CONTACT_MODE:-closest}"          # back to ClosestPoint

# --- previously-good per-shape bounds (your confirmed values) ---
export DISP_sphere="${DISP_sphere:-34.8}"
export STEPS_sphere="${STEPS_sphere:-348}"

export DISP_oblate="${DISP_oblate:-19.7}"
export STEPS_oblate="${STEPS_oblate:-197}"

export DISP_cone="${DISP_cone:-15.9}"
export STEPS_cone="${STEPS_cone:-159}"

export DISP_capcone="${DISP_capcone:-17.7}"
export STEPS_capcone="${STEPS_capcone:-177}"

export DISP_doublecone="${DISP_doublecone:-19.3}"
export STEPS_doublecone="${STEPS_doublecone:-193}"

# Skip prolate and cap — still being diagnosed.
export SHAPES="${SHAPES:-sphere oblate cone capcone doublecone}"

# Don't double-scale displacements (overrides above are already final).
export DISP_SCALE="${DISP_SCALE:-1.0}"

echo "==== Revert sweep settings ===="
echo "  CONTACT_MODE       = $CONTACT_MODE"
echo "  DISP_SCALE         = $DISP_SCALE"
echo "  SHAPES             = $SHAPES"
echo "  per-shape DISP / STEPS:"
for s in $SHAPES; do
    d_var="DISP_$s"; n_var="STEPS_$s"
    echo "    $s  DISP=${!d_var}  STEPS=${!n_var}"
done
echo

./run_final_sweep.sh
