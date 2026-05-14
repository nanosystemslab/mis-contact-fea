#!/usr/bin/env bash
# Run a sequence of YAML configs through the full pipeline.
#
# Usage:
#   ./scripts/run_configs.sh examples/01_seven_builtins/*_push.yaml
#   ./scripts/run_configs.sh examples/01_seven_builtins/sphere_push.yaml \
#                            examples/01_seven_builtins/cap_push.yaml
#
# Each config runs through run_config.sh; failures are logged but don't
# abort the sweep, so one shape stalling won't kill the rest.
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

if [ "$#" -lt 1 ]; then
    echo "usage: $0 path/to/config1.yaml [config2.yaml ...]" >&2
    exit 1
fi

CONFIGS=("$@")
SUCCEEDED=()
FAILED=()

echo "==== Running ${#CONFIGS[@]} configs ===="
for cfg in "${CONFIGS[@]}"; do
    echo
    echo "########################################################################"
    echo "##  $cfg"
    echo "########################################################################"
    if ./scripts/run_config.sh "$cfg"; then
        SUCCEEDED+=("$cfg")
    else
        FAILED+=("$cfg")
        echo "!! $cfg FAILED, continuing"
    fi
done

echo
echo "==== Sweep complete ===="
echo "  succeeded (${#SUCCEEDED[@]}):"
for c in "${SUCCEEDED[@]}"; do echo "    $c"; done
if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "  failed (${#FAILED[@]}):"
    for c in "${FAILED[@]}"; do echo "    $c"; done
fi
