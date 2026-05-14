#!/usr/bin/env bash
# Launch the 2D alignment GUI in a browser tab.
#
# Requires streamlit (install via `pip install -e ".[gui]"` from project root).
set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

if ! command -v streamlit >/dev/null 2>&1; then
    echo "Streamlit not found. Install with:" >&2
    echo "    pip install -e \".[gui]\"" >&2
    echo "  or" >&2
    echo "    pip install streamlit" >&2
    exit 1
fi

PORT="${PORT:-8501}"
streamlit run src/mis_contact_fea/gui/align_2d.py \
    --server.port "$PORT" \
    --browser.gatherUsageStats false \
    "$@"
