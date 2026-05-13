#!/usr/bin/env bash
# Build the dolfinx + dolfinx_contact (asimov-contact) runtime image.
#
# The image is tagged BOTH with the local name used by the pipeline
# scripts (`dolfinx-contact:local`) AND with the Docker Hub coordinates
# (`mattnakamura/dolfinx-contact:v0.9.0` and `:latest`), so a single
# build serves both local use and a subsequent `docker push`.
#
# Usage:
#   ./docker/build_asimov_image.sh                 # default tags
#   ./docker/build_asimov_image.sh <base_image> <local_tag> <hub_repo> <hub_tag>
#
# To publish to Docker Hub after building:
#   docker login -u mattnakamura
#   docker push mattnakamura/dolfinx-contact:v0.9.0
#   docker push mattnakamura/dolfinx-contact:latest
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

BASE_IMAGE="${1:-ghcr.io/fenics/dolfinx/dolfinx:v0.9.0}"
LOCAL_TAG="${2:-dolfinx-contact:local}"
HUB_REPO="${3:-mattnakamura/dolfinx-contact}"
HUB_TAG="${4:-v0.9.0}"

echo "=== Building dolfinx-contact image ==="
echo "  base:       $BASE_IMAGE"
echo "  local tag:  $LOCAL_TAG"
echo "  hub tag:    $HUB_REPO:$HUB_TAG (and :latest)"
echo

docker build \
  --build-arg BASE_IMAGE="$BASE_IMAGE" \
  -f docker/Dockerfile.asimov \
  -t "$LOCAL_TAG" \
  -t "$HUB_REPO:$HUB_TAG" \
  -t "$HUB_REPO:latest" \
  .

echo
echo "Done. To publish to Docker Hub:"
echo "  docker login -u ${HUB_REPO%%/*}"
echo "  docker push $HUB_REPO:$HUB_TAG"
echo "  docker push $HUB_REPO:latest"
