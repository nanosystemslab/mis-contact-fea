#!/usr/bin/env bash
# Rebuild the comparative force / grid / R-vs-I plots from whatever is
# currently sitting in runs_all_shapes/ (push) and runs_retention/.
#
# Use this AFTER running several single-shape simulations to assemble
# the multi-shape comparison without re-simulating anything.
#
# Override roots via env:
#   PUSH_ROOT=runs_smoke_test ./scripts/build_comparison_plots.sh
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

PUSH_ROOT="${PUSH_ROOT:-runs_all_shapes}"
RETENTION_ROOT="${RETENTION_ROOT:-runs_retention}"

generated=()

if [ -d "$PUSH_ROOT" ]; then
    echo "=== Push comparison ($PUSH_ROOT) ==="
    python3 -m mis_contact_fea.postproc.compare_forces --abs --mode contact-aligned \
        --contact-threshold 1 --root "$PUSH_ROOT" \
        --out force_comparison_aligned.png && generated+=("force_comparison_aligned.png") || true
    python3 -m mis_contact_fea.postproc.compare_shapes --root "$PUSH_ROOT" \
        --out comparison_grid.png && generated+=("comparison_grid.png") || true
    echo
fi

if [ -d "$RETENTION_ROOT" ]; then
    echo "=== Retention comparison ($RETENTION_ROOT) ==="
    python3 -m mis_contact_fea.postproc.compare_forces --abs --mode contact-aligned \
        --contact-threshold 1 --root "$RETENTION_ROOT" \
        --out force_comparison_retention.png && generated+=("force_comparison_retention.png") || true
    python3 -m mis_contact_fea.postproc.compare_shapes --root "$RETENTION_ROOT" \
        --out comparison_grid_retention.png && generated+=("comparison_grid_retention.png") || true
    echo
fi

if [ -d "$PUSH_ROOT" ] && [ -d "$RETENTION_ROOT" ]; then
    echo "=== Smoothed + R/I analysis ==="
    python3 -m mis_contact_fea.postproc.analyze && {
        generated+=("force_smoothed_overlay.png" "force_smoothed_push.png"
                    "force_smoothed_retention.png" "force_ratio.png" "force_peaks.csv")
    } || true
    echo
fi

if [ ${#generated[@]} -gt 0 ]; then
    echo "Generated:"
    for f in "${generated[@]}"; do echo "  $f"; done
    # Open only the PNGs (skip CSV)
    open_targets=()
    for f in "${generated[@]}"; do
        [[ "$f" == *.png ]] && open_targets+=("$f")
    done
    [ ${#open_targets[@]} -gt 0 ] && open "${open_targets[@]}" 2>/dev/null || true
else
    echo "Nothing generated — neither $PUSH_ROOT/ nor $RETENTION_ROOT/ exists."
fi
