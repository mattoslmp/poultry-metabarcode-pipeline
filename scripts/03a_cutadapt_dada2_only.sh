#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${QIIME_DIR}" "${LOG_DIR}" "${PROJECT_DIR}/tmp_qiime"
exec > >(tee -a "${LOG_DIR}/03a_cutadapt_dada2_only.log") 2>&1

export TMPDIR="${TMPDIR:-${PROJECT_DIR}/tmp_qiime}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

if [[ -s "${MANUAL_METADATA_TSV}" ]] && [[ $(grep -vc '^#' "${MANUAL_METADATA_TSV}") -gt 1 ]]; then
  METADATA_TSV="${MANUAL_METADATA_TSV}"
else
  METADATA_TSV="${AUTO_METADATA_TSV}"
fi

[[ -s "${QIIME_DIR}/demux-prinseq.qza" ]] || { echo "Missing ${QIIME_DIR}/demux-prinseq.qza. Run step 02 first."; exit 1; }
[[ -s "${METADATA_TSV}" ]] || { echo "Missing metadata: ${METADATA_TSV}"; exit 1; }

echo "[INFO] Removing 515F/806R primers and reverse complements with q2-cutadapt..."
qiime cutadapt trim-single \
  --i-demultiplexed-sequences "${QIIME_DIR}/demux-prinseq.qza" \
  --p-front "^${PRIMER_515F}" \
  --p-front "^${PRIMER_806R}" \
  --p-front "^${PRIMER_515F_RC}" \
  --p-front "^${PRIMER_806R_RC}" \
  --p-overlap "${CUTADAPT_OVERLAP}" \
  --p-error-rate "${CUTADAPT_ERROR_RATE}" \
  --p-match-adapter-wildcards true \
  --p-match-read-wildcards true \
  --p-cores "${THREADS}" \
  --o-trimmed-sequences "${QIIME_DIR}/demux-prinseq-cutadapt.qza" \
  --verbose

echo "[INFO] Denoising with DADA2 single-end..."
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${QIIME_DIR}/demux-prinseq-cutadapt.qza" \
  --p-trim-left "${DADA2_TRIM_LEFT}" \
  --p-trunc-len "${DADA2_TRUNC_LEN}" \
  --p-chimera-method "${DADA2_CHIMERA_METHOD}" \
  --p-n-threads "${THREADS}" \
  --o-table "${QIIME_DIR}/table.qza" \
  --o-representative-sequences "${QIIME_DIR}/rep-seqs.qza" \
  --o-denoising-stats "${QIIME_DIR}/denoising-stats.qza" \
  --verbose

qiime feature-table summarize \
  --i-table "${QIIME_DIR}/table.qza" \
  --m-sample-metadata-file "${METADATA_TSV}" \
  --o-visualization "${QIIME_DIR}/table.qzv"

qiime metadata tabulate \
  --m-input-file "${QIIME_DIR}/denoising-stats.qza" \
  --o-visualization "${QIIME_DIR}/denoising-stats.qzv"

echo "[DONE] Cutadapt and DADA2 finished. Continue with scripts/03b_resume_after_dada2_lowmem_taxonomy.sh"
