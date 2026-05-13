#!/bin/bash
# Build DOLFINx container with dolfinx_contact support

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_CONTAINER="${BASE_CONTAINER:-${HOME}/Optimization_Framework/containers/dolfinx_v0.9.0.sif}"
OUTPUT_CONTAINER="${OUTPUT_CONTAINER:-${HOME}/mis-contact-fea/containers/dolfinx_v0.9.0_contact.sif}"
DEF_FILE="${SCRIPT_DIR}/dolfinx_v0.9.0_contact.def"

echo "=== Building DOLFINx Container with dolfinx_contact ==="
echo "Base container:   ${BASE_CONTAINER}"
echo "Output container: ${OUTPUT_CONTAINER}"
echo "Definition file:  ${DEF_FILE}"
echo ""

if [ ! -f "${BASE_CONTAINER}" ]; then
    echo "ERROR: Base container not found: ${BASE_CONTAINER}"
    echo ""
    echo "Set BASE_CONTAINER or ensure the base DOLFINx container exists."
    exit 1
fi

if [ ! -f "${DEF_FILE}" ]; then
    echo "ERROR: Definition file not found: ${DEF_FILE}"
    exit 1
fi

if [ -f "${OUTPUT_CONTAINER}" ]; then
    echo "WARNING: Output container already exists: ${OUTPUT_CONTAINER}"
    read -p "Overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
    rm -f "${OUTPUT_CONTAINER}"
fi

if singularity help build | grep -q "\-\-fakeroot"; then
    FAKEROOT_FLAG="--fakeroot"
else
    FAKEROOT_FLAG=""
fi

singularity build ${FAKEROOT_FLAG} "${OUTPUT_CONTAINER}" "${DEF_FILE}"

echo ""
echo "=== Build Complete ==="
echo "Container: ${OUTPUT_CONTAINER}"
echo ""

echo "=== Testing dolfinx_contact ==="
singularity exec "${OUTPUT_CONTAINER}" python3 << 'EOF'
import sys
print(f"Python: {sys.version}")
print()

try:
    import dolfinx_contact
    print("✓ dolfinx_contact available")
except ImportError as e:
    print(f"✗ dolfinx_contact not available: {e}")
    sys.exit(1)

try:
    import dolfinx
    print(f"✓ DOLFINx {dolfinx.__version__} available")
except ImportError as e:
    print(f"✗ DOLFINx not available: {e}")
    sys.exit(1)

print()
print("Container is ready to use!")
EOF

echo ""
echo "=== Usage ==="
echo "Set CONTAINER to use the new image:"
echo "  CONTAINER=\"${OUTPUT_CONTAINER}\""
