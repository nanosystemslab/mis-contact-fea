#!/bin/bash
#SBATCH --job-name=contact_fea
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=mtdsn@hawaii.edu
#SBATCH --output=logs/contact_%j.log
#SBATCH --error=logs/contact_%j.err

set -e
set -u

echo "======================================================================"
echo "Contact FEA - Job ${SLURM_JOB_ID}"
echo "======================================================================"
echo "Started at: $(date)"
echo "Running on node: $(hostname)"
echo ""

PROJECT_DIR="${PROJECT_DIR:-${HOME}/basic_fea}"
CONTAINER="${CONTAINER:-${HOME}/basic_fea/containers/dolfinx_v0.9.0_contact.sif}"
OUTPUT_DIR="${OUTPUT_DIR:-${HOME}/scratch/contact_${SLURM_JOB_ID}}"

BASE_STEP="${BASE_STEP:-models/N1004_base.STEP}"
STRUCTURE_STEP="${STRUCTURE_STEP:-models/N1004_structure.STEP}"

CLEARANCE_MM="${CLEARANCE_MM:-1.0}"
MESH_SIZE_MM="${MESH_SIZE_MM:-2.0}"
SURFACE_TOL_MM="${SURFACE_TOL_MM:-2.0}"
DISP_MM="${DISP_MM:-auto}"
STEPS="${STEPS:-20}"
FRIC="${FRIC:-0.0}"
GAMMA_SCALE="${GAMMA_SCALE:-1.0}"
RELAX="${RELAX:-0.5}"
NEWTON_MAX_IT="${NEWTON_MAX_IT:-80}"
NEWTON_RTOL="${NEWTON_RTOL:-1e-6}"
NEWTON_ATOL="${NEWTON_ATOL:-1e-7}"
KSP_RTOL="${KSP_RTOL:-1e-8}"
KSP_ATOL="${KSP_ATOL:-1e-12}"
KSP_MAX_IT="${KSP_MAX_IT:-2000}"
KSP_TYPE="${KSP_TYPE:-fgmres}"
SKIP_MESH="${SKIP_MESH:-0}"
MESH_DIR="${MESH_DIR:-}"
PC_TYPE="${PC_TYPE:-gamg}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

module purge
module load tools/Singularity/latest 2>/dev/null || \
module load tools/Apptainer/latest 2>/dev/null

if [ ! -f "${CONTAINER}" ]; then
    echo "ERROR: Container not found at ${CONTAINER}"
    exit 1
fi

CONTAINER_CMD=$(command -v singularity || command -v apptainer)
RANKS="${RANKS:-${SLURM_NTASKS:-1}}"
export OMP_NUM_THREADS=1

echo "Container: ${CONTAINER}"
echo "Project:   ${PROJECT_DIR}"
echo "Output:    ${OUTPUT_DIR}"
echo "Skip mesh: ${SKIP_MESH}"
echo "Mesh dir:  ${MESH_DIR:-${OUTPUT_DIR}/mesh_out}"
echo "Base:      ${BASE_STEP}"
echo "Struct:    ${STRUCTURE_STEP}"
echo "Clearance: ${CLEARANCE_MM} mm"
echo "Surf tol:  ${SURFACE_TOL_MM} mm"
echo "Disp:      ${DISP_MM} mm"
echo "Steps:     ${STEPS}"
echo "Friction:  ${FRIC}"
echo "Gamma:     ${GAMMA_SCALE}"
echo "Relax:     ${RELAX}"
echo "Newton:    max_it=${NEWTON_MAX_IT} rtol=${NEWTON_RTOL} atol=${NEWTON_ATOL}"
echo "KSP:       type=${KSP_TYPE} rtol=${KSP_RTOL} atol=${KSP_ATOL} max_it=${KSP_MAX_IT}"
echo "PC:        ${PC_TYPE}"
echo ""

export XDG_CACHE_HOME="${OUTPUT_DIR}/.xdg_cache_fenicsx"
mkdir -p "${XDG_CACHE_HOME}/fenics"
export DOLFINX_JIT_TIMEOUT="${DOLFINX_JIT_TIMEOUT:-3600}"

# Determine mesh directory (use existing or create new)
if [ -n "${MESH_DIR}" ]; then
    ACTUAL_MESH_DIR="${MESH_DIR}"
else
    ACTUAL_MESH_DIR="${OUTPUT_DIR}/mesh_out"
fi

${CONTAINER_CMD} exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${OUTPUT_DIR}:${OUTPUT_DIR}" \
    ${MESH_DIR:+--bind "${MESH_DIR}:${MESH_DIR}"} \
    --pwd "${PROJECT_DIR}" \
    "${CONTAINER}" \
    bash -c "
export PYTHONUNBUFFERED=1

if [ '${SKIP_MESH}' != '1' ]; then
  python3 -u src/build_contact_mesh.py \
    --base-step ${BASE_STEP} \
    --structure-step ${STRUCTURE_STEP} \
    --out-dir ${ACTUAL_MESH_DIR} \
    --clearance-mm ${CLEARANCE_MM} \
    --mesh-size-mm ${MESH_SIZE_MM} \
    --surface-tol-mm ${SURFACE_TOL_MM} \
    --flip-xy
else
  echo 'Skipping mesh generation (SKIP_MESH=1), using existing mesh at ${ACTUAL_MESH_DIR}'
fi

echo 'Starting contact solve...'
DISP_ARG=''
if [ '${DISP_MM}' != 'auto' ]; then
  DISP_ARG=\"--disp-mm ${DISP_MM}\"
fi

if [ '${RANKS}' -le 1 ]; then
  python3 -u src/run_contact_sim.py \
    --mesh-dir ${ACTUAL_MESH_DIR} \
    --out-dir ${OUTPUT_DIR}/results_contact \
    --steps ${STEPS} \
    \${DISP_ARG} \
    --fric ${FRIC} \
    --gamma-scale ${GAMMA_SCALE} \
    --relax ${RELAX} \
    --newton-max-it ${NEWTON_MAX_IT} \
    --newton-rtol ${NEWTON_RTOL} \
    --newton-atol ${NEWTON_ATOL} \
    --ksp-type ${KSP_TYPE} \
    --pc-type ${PC_TYPE} \
    --ksp-rtol ${KSP_RTOL} \
    --ksp-atol ${KSP_ATOL} \
    --ksp-max-it ${KSP_MAX_IT}
else
  env HYDRA_LAUNCHER=fork OMP_NUM_THREADS=1 mpiexec -n ${RANKS} python3 -u src/run_contact_sim.py \
    --mesh-dir ${ACTUAL_MESH_DIR} \
    --out-dir ${OUTPUT_DIR}/results_contact \
    --steps ${STEPS} \
    \${DISP_ARG} \
    --fric ${FRIC} \
    --gamma-scale ${GAMMA_SCALE} \
    --relax ${RELAX} \
    --newton-max-it ${NEWTON_MAX_IT} \
    --newton-rtol ${NEWTON_RTOL} \
    --newton-atol ${NEWTON_ATOL} \
    --ksp-type ${KSP_TYPE} \
    --pc-type ${PC_TYPE} \
    --ksp-rtol ${KSP_RTOL} \
    --ksp-atol ${KSP_ATOL} \
    --ksp-max-it ${KSP_MAX_IT}
fi
"

EXIT_CODE=$?
echo ""
echo "======================================================================"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Contact FEA finished successfully"
    echo "Results: ${OUTPUT_DIR}/results_contact"
else
    echo "✗ Contact FEA failed (exit code: ${EXIT_CODE})"
fi
echo "======================================================================"

exit ${EXIT_CODE}
