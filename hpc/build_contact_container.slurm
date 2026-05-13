#!/bin/bash
#SBATCH --job-name=build_contact_container
#SBATCH --partition=sandbox
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=mtdsn@hawaii.edu
#SBATCH --output=logs/build_contact_container_%j.log
#SBATCH --error=logs/build_contact_container_%j.err

set -e
set -u

echo "======================================================================"
echo "Build dolfinx_contact Container - Job ${SLURM_JOB_ID}"
echo "======================================================================"
echo "Started at: $(date)"
echo "Running on node: $(hostname)"
echo ""

module purge
module load tools/Singularity/latest 2>/dev/null || \
module load tools/Apptainer/latest 2>/dev/null || {
    echo "ERROR: Neither Singularity nor Apptainer available"
    exit 1
}

if command -v singularity &> /dev/null; then
    RUNTIME="singularity"
elif command -v apptainer &> /dev/null; then
    RUNTIME="apptainer"
else
    echo "ERROR: Neither 'singularity' nor 'apptainer' found in PATH"
    exit 1
fi

PROJECT_DIR="${PROJECT_DIR:-${HOME}/basic_fea}"
BASE_CONTAINER="${BASE_CONTAINER:-${HOME}/Optimization_Framework/containers/dolfinx_v0.9.0.sif}"
OUTPUT_CONTAINER="${OUTPUT_CONTAINER:-${HOME}/basic_fea/containers/dolfinx_v0.9.0_contact.sif}"

mkdir -p "${PROJECT_DIR}/logs"

echo "Runtime:   ${RUNTIME}"
echo "Project:   ${PROJECT_DIR}"
echo "Base SIF:  ${BASE_CONTAINER}"
echo "Output:    ${OUTPUT_CONTAINER}"
echo ""

if [ ! -f "${BASE_CONTAINER}" ]; then
    echo "ERROR: Base container not found at ${BASE_CONTAINER}"
    exit 1
fi

cd "${PROJECT_DIR}/containers"

export BASE_CONTAINER
export OUTPUT_CONTAINER

echo "Building container..."
${RUNTIME} exec "${BASE_CONTAINER}" /bin/true || true

bash build_contact_container.sh

echo ""
echo "======================================================================"
echo "✓ Container build finished"
echo "======================================================================"
