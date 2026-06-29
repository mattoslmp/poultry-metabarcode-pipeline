#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${QIIME_DIR}" "${EXPORTED_DIR}" "${LOG_DIR}" "${PROJECT_DIR}/tmp_qiime"
exec > >(tee -a "${LOG_DIR}/03c_resume_after_taxonomy_no_ordination.log") 2>&1

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
[[ -s "${QIIME_DIR}/taxonomy.qza" ]] || { echo "Missing ${QIIME_DIR}/taxonomy.qza"; exit 1; }
[[ -s "${QIIME_DIR}/denoising-stats.qza" ]] || { echo "Missing ${QIIME_DIR}/denoising-stats.qza"; exit 1; }
[[ -s "${METADATA_TSV}" ]] || { echo "Missing metadata: ${METADATA_TSV}"; exit 1; }

echo "[INFO] Resuming after taxonomy. This script avoids core-metrics beta ordination."
echo "[INFO] Rarefaction depth: ${RAREFY_DEPTH}"
echo "[INFO] Metadata: ${METADATA_TSV}"

for level in 3 6; do
  if [[ ! -s "${QIIME_DIR}/taxa_level${level}.qza" ]]; then
    qiime taxa collapse \
      --i-table "${QIIME_DIR}/table.qza" \
      --i-taxonomy "${QIIME_DIR}/taxonomy.qza" \
      --p-level "${level}" \
      --o-collapsed-table "${QIIME_DIR}/taxa_level${level}.qza"
  else
    echo "[SKIP] ${QIIME_DIR}/taxa_level${level}.qza exists"
  fi
done

if [[ ! -s "${QIIME_DIR}/taxa-bar-plots.qzv" ]]; then
  qiime taxa barplot \
    --i-table "${QIIME_DIR}/table.qza" \
    --i-taxonomy "${QIIME_DIR}/taxonomy.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "${QIIME_DIR}/taxa-bar-plots.qzv"
fi

qiime feature-table rarefy \
  --i-table "${QIIME_DIR}/table.qza" \
  --p-sampling-depth "${RAREFY_DEPTH}" \
  --o-rarefied-table "${QIIME_DIR}/rarefied_table_${RAREFY_DEPTH}.qza"

qiime diversity alpha \
  --i-table "${QIIME_DIR}/rarefied_table_${RAREFY_DEPTH}.qza" \
  --p-metric shannon \
  --o-alpha-diversity "${QIIME_DIR}/shannon_vector.qza"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "${QIIME_DIR}/shannon_vector.qza" \
  --m-metadata-file "${METADATA_TSV}" \
  --o-visualization "${QIIME_DIR}/shannon_group_significance.qzv" || true

rm -rf \
  "${EXPORTED_DIR}/feature_table" \
  "${EXPORTED_DIR}/rarefied_table" \
  "${EXPORTED_DIR}/denoising_stats" \
  "${EXPORTED_DIR}/taxa_level3" \
  "${EXPORTED_DIR}/taxa_level6" \
  "${EXPORTED_DIR}/taxonomy" \
  "${EXPORTED_DIR}/shannon"

qiime tools export --input-path "${QIIME_DIR}/table.qza" --output-path "${EXPORTED_DIR}/feature_table"
qiime tools export --input-path "${QIIME_DIR}/rarefied_table_${RAREFY_DEPTH}.qza" --output-path "${EXPORTED_DIR}/rarefied_table"
qiime tools export --input-path "${QIIME_DIR}/denoising-stats.qza" --output-path "${EXPORTED_DIR}/denoising_stats"
qiime tools export --input-path "${QIIME_DIR}/taxa_level3.qza" --output-path "${EXPORTED_DIR}/taxa_level3"
qiime tools export --input-path "${QIIME_DIR}/taxa_level6.qza" --output-path "${EXPORTED_DIR}/taxa_level6"
qiime tools export --input-path "${QIIME_DIR}/taxonomy.qza" --output-path "${EXPORTED_DIR}/taxonomy"
qiime tools export --input-path "${QIIME_DIR}/shannon_vector.qza" --output-path "${EXPORTED_DIR}/shannon"

biom convert -i "${EXPORTED_DIR}/feature_table/feature-table.biom" -o "${EXPORTED_DIR}/feature_table.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/rarefied_table/feature-table.biom" -o "${EXPORTED_DIR}/rarefied_table.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/taxa_level3/feature-table.biom" -o "${EXPORTED_DIR}/taxa_level3.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/taxa_level6/feature-table.biom" -o "${EXPORTED_DIR}/taxa_level6.tsv" --to-tsv
cp "${METADATA_TSV}" "${EXPORTED_DIR}/sample-metadata.used.tsv"
echo "${RAREFY_DEPTH}" > "${EXPORTED_DIR}/rarefaction_depth.txt"

echo "[DONE] Resume after taxonomy complete without beta ordinations."
