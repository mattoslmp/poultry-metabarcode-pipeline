#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_DIR}/config/pipeline_config.sh"

TAXONOMY_JOBS="${TAXONOMY_JOBS:-1}"
READS_PER_BATCH="${READS_PER_BATCH:-50}"

mkdir -p "${QIIME_DIR}" "${EXPORTED_DIR}" "${LOG_DIR}" "${PROJECT_DIR}/tmp_qiime"
exec > >(tee -a "${LOG_DIR}/03b_resume_lowmem_taxonomy.log") 2>&1

export TMPDIR="${TMPDIR:-${PROJECT_DIR}/tmp_qiime}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -s "${MANUAL_METADATA_TSV}" ]] && [[ $(grep -vc '^#' "${MANUAL_METADATA_TSV}") -gt 1 ]]; then
  METADATA_TSV="${MANUAL_METADATA_TSV}"
else
  METADATA_TSV="${AUTO_METADATA_TSV}"
fi

[[ -s "${QIIME_DIR}/table.qza" ]] || { echo "Missing ${QIIME_DIR}/table.qza"; exit 1; }
[[ -s "${QIIME_DIR}/rep-seqs.qza" ]] || { echo "Missing ${QIIME_DIR}/rep-seqs.qza"; exit 1; }
[[ -s "${QIIME_DIR}/denoising-stats.qza" ]] || { echo "Missing ${QIIME_DIR}/denoising-stats.qza"; exit 1; }
[[ -s "${CLASSIFIER_QZA}" ]] || { echo "Missing classifier: ${CLASSIFIER_QZA}"; exit 1; }
[[ -s "${METADATA_TSV}" ]] || { echo "Missing metadata: ${METADATA_TSV}"; exit 1; }

echo "[INFO] TMPDIR=$TMPDIR"
echo "[INFO] TAXONOMY_JOBS=$TAXONOMY_JOBS"
echo "[INFO] READS_PER_BATCH=$READS_PER_BATCH"
echo "[INFO] Running low-memory taxonomy classification only."

qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER_QZA}" \
  --i-reads "${QIIME_DIR}/rep-seqs.qza" \
  --p-n-jobs "${TAXONOMY_JOBS}" \
  --p-reads-per-batch "${READS_PER_BATCH}" \
  --o-classification "${QIIME_DIR}/taxonomy.qza" \
  --verbose

for level in 3 6; do
  qiime taxa collapse \
    --i-table "${QIIME_DIR}/table.qza" \
    --i-taxonomy "${QIIME_DIR}/taxonomy.qza" \
    --p-level "${level}" \
    --o-collapsed-table "${QIIME_DIR}/taxa_level${level}.qza"
done

qiime taxa barplot \
  --i-table "${QIIME_DIR}/table.qza" \
  --i-taxonomy "${QIIME_DIR}/taxonomy.qza" \
  --m-metadata-file "${METADATA_TSV}" \
  --o-visualization "${QIIME_DIR}/taxa-bar-plots.qzv"

echo "[DONE] Low-memory taxonomy classification and L3/L6 taxonomic tables complete. Continue with scripts/03c_resume_after_taxonomy_no_ordination.sh"
