#!/bin/bash
#SBATCH --job-name=infercnv
#SBATCH --account=project_xxxxxxx
#SBATCH --partition=small
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-1%4
#SBATCH --output=cnv_inference/logs/%x_%A_%a.out
#SBATCH --error=cnv_inference/logs/%x_%A_%a.err

set -euo pipefail

module load r-env

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

MANIFEST="${REPO_ROOT}/copy-number-variations/external-cohort/cnv_inference/infercnv_input/tumor-stroma/manifest.tsv"
R_SCRIPT="${REPO_ROOT}/copy-number-variations/scripts/run_infercnv_sample.R"
GENE_ORDER="${REPO_ROOT}/copy-number-variations/files/gene_order_file.tsv"

mkdir -p logs
mkdir -p "$(dirname "$MANIFEST")"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")

IFS=$'\t' read -r SAMPLE COUNTS_H5 ANNOTATIONS_TSV OUT_DIR CUTOFF HMM REF_GROUPS <<< "$LINE"

mkdir -p "$OUT_DIR"

if [[ -n "${PROJECT_R_LIB:-}" ]]; then
  export R_LIBS_USER="${PROJECT_R_LIB}:${R_LIBS_USER:-}"
fi

echo "Running sample: $SAMPLE"
echo "Counts: $COUNTS_H5"
echo "Annotations: $ANNOTATIONS_TSV"
echo "Out dir: $OUT_DIR"
echo "Reference groups: $REF_GROUPS"
echo "HMM: $HMM"

cmd=(
  Rscript "$R_SCRIPT"
  --counts_h5 "$COUNTS_H5"
  --annotations_tsv "$ANNOTATIONS_TSV"
  --gene_order_file "$GENE_ORDER"
  --out_dir "$OUT_DIR"
  --cutoff "$CUTOFF"
  --cluster_by_groups TRUE
  --denoise TRUE
  --num_threads 8
)

# ---------- optional ref_groups ----------
if [[ -n "${REF_GROUPS:-}" ]]; then
  cmd+=("--ref_groups=$REF_GROUPS")
fi

# ---------- optional HMM ----------
if [[ -n "${HMM:-}" ]]; then
  case "$HMM" in
    TRUE|FALSE)
      cmd+=(--HMM "$HMM")
      ;;
    *)
      echo "Error: HMM must be TRUE, FALSE, or empty in manifest, got: '$HMM'" >&2
      exit 1
      ;;
  esac
fi

printf 'Command arg: <%s>\n' "${cmd[@]}"

"${cmd[@]}"
