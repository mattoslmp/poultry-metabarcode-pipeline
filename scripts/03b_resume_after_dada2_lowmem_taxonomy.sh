#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_DIR}/config/pipeline_config.sh"

TAXONOMY_JOBS="${TAXONOMY_JOBS:-1}"
READS_PER_BATCH="${READS_PER_BATCH:-100}"

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
echo "[INFO] Running low-memory taxonomy classification..."

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

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME_DIR}/rep-seqs.qza" \
  --p-n-threads 1 \
  --o-alignment "${QIIME_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME_DIR}/rooted-tree.qza"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${QIIME_DIR}/rooted-tree.qza" \
  --i-table "${QIIME_DIR}/table.qza" \
  --p-sampling-depth "${RAREFY_DEPTH}" \
  --m-metadata-file "${METADATA_TSV}" \
  --output-dir "${QIIME_DIR}/core-metrics-phylogenetic-${RAREFY_DEPTH}"

qiime diversity alpha-rarefaction \
  --i-table "${QIIME_DIR}/table.qza" \
  --i-phylogeny "${QIIME_DIR}/rooted-tree.qza" \
  --p-max-depth "${RAREFY_DEPTH}" \
  --m-metadata-file "${METADATA_TSV}" \
  --o-visualization "${QIIME_DIR}/alpha_rarefaction_max_${RAREFY_DEPTH}.qzv"

rm -rf \
  "${EXPORTED_DIR}/feature_table" \
  "${EXPORTED_DIR}/rarefied_table" \
  "${EXPORTED_DIR}/denoising_stats" \
  "${EXPORTED_DIR}/taxa_level3" \
  "${EXPORTED_DIR}/taxa_level6" \
  "${EXPORTED_DIR}/taxonomy" \
  "${EXPORTED_DIR}/shannon"

qiime tools export --input-path "${QIIME_DIR}/table.qza" --output-path "${EXPORTED_DIR}/feature_table"
qiime tools export --input-path "${QIIME_DIR}/core-metrics-phylogenetic-${RAREFY_DEPTH}/rarefied_table.qza" --output-path "${EXPORTED_DIR}/rarefied_table"
qiime tools export --input-path "${QIIME_DIR}/denoising-stats.qza" --output-path "${EXPORTED_DIR}/denoising_stats"
qiime tools export --input-path "${QIIME_DIR}/taxa_level3.qza" --output-path "${EXPORTED_DIR}/taxa_level3"
qiime tools export --input-path "${QIIME_DIR}/taxa_level6.qza" --output-path "${EXPORTED_DIR}/taxa_level6"
qiime tools export --input-path "${QIIME_DIR}/taxonomy.qza" --output-path "${EXPORTED_DIR}/taxonomy"
qiime tools export --input-path "${QIIME_DIR}/core-metrics-phylogenetic-${RAREFY_DEPTH}/shannon_vector.qza" --output-path "${EXPORTED_DIR}/shannon"

biom convert -i "${EXPORTED_DIR}/feature_table/feature-table.biom" -o "${EXPORTED_DIR}/feature_table.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/rarefied_table/feature-table.biom" -o "${EXPORTED_DIR}/rarefied_table.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/taxa_level3/feature-table.biom" -o "${EXPORTED_DIR}/taxa_level3.tsv" --to-tsv
biom convert -i "${EXPORTED_DIR}/taxa_level6/feature-table.biom" -o "${EXPORTED_DIR}/taxa_level6.tsv" --to-tsv
cp "${METADATA_TSV}" "${EXPORTED_DIR}/sample-metadata.used.tsv"
echo "${RAREFY_DEPTH}" > "${EXPORTED_DIR}/rarefaction_depth.txt"

echo "[DONE] Low-memory taxonomy resume complete."
