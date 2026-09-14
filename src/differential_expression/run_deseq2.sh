#!/bin/bash
#SBATCH --job-name=deseq2_visium
#SBATCH --account=project_xxxxxxx
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --output=logs/deseq2_%j.out
#SBATCH --error=logs/deseq2_%j.err
#SBATCH --mail-type=END,FAIL

# ── Environment ───────────────────────────────────────────────────────────────
# Edit these paths for your system. ENV_PATH should point at a Python
# environment built from the repo's pyproject.toml / uv.lock.
module load python-data

if [[ -n "${ENV_PATH:-}" ]]; then
  export PATH="${ENV_PATH}:${PATH}"
fi


REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ANNDATA_PATH="${REPO_ROOT}/quality_control/primary-cohort/adata.h5ad"
OUT_DEGS="${REPO_ROOT}/differential_expression/primary-cohort/deseq2_results_batch"
SCRIPT_DIR="${REPO_ROOT}/src/differential_expression"

# ── Create directories ────────────────────────────────────────────────────────
mkdir -p logs
mkdir -p "${OUT_DEGS}"

# ── Print job info to log ─────────────────────────────────────────────────────
echo "Job ID        : ${SLURM_JOB_ID}"
echo "Node          : ${SLURMD_NODENAME}"
echo "CPUs          : ${SLURM_CPUS_PER_TASK}"
echo "Memory        : ${SLURM_MEM_PER_NODE} MB"
echo "Environment:  : ${ENV_PATH:-<system>}"
echo "AnnData       : ${ANNDATA_PATH}"
echo "Output dir    : ${OUT_DEGS}"
echo "Start time    : $(date)"
echo "─────────────────────────────────────────────────────"

# ── Run ───────────────────────────────────────────────────────────────────────
python "${SCRIPT_DIR}/run_deseq2.py" \
    --anndata "${ANNDATA_PATH}" \
    --outdir "${OUT_DEGS}" \
    --min-count 10 \
    --min-samples 2 \
    --padj 0.05 \
    --lfc 1.0 \
    --ncpus "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "─────────────────────────────────────────────────────"
echo "End time      : $(date)"
echo "Exit code     : ${EXIT_CODE}"
exit ${EXIT_CODE}
