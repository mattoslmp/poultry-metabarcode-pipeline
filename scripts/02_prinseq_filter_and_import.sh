#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${FILTERED_DIR}" "${METADATA_DIR}" "${QIIME_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/02_prinseq_filter_and_import.log") 2>&1

if command -v prinseq-lite.pl >/dev/null 2>&1; then
  PRINSEQ_CMD="prinseq-lite.pl"
elif command -v prinseq-lite >/dev/null 2>&1; then
  PRINSEQ_CMD="prinseq-lite"
else
  echo "[ERROR] PRINSEQ-lite was not found in PATH."
  exit 1
fi

echo "[INFO] Using PRINSEQ command: ${PRINSEQ_CMD}"
echo "[INFO] Filtering: min_len=${PRINSEQ_MIN_LEN}; min_qual_mean=${PRINSEQ_MIN_QUAL_MEAN}"

shopt -s nullglob
fastqs=("${RAW_DIR}"/*.fastq.gz "${RAW_DIR}"/*.fq.gz "${RAW_DIR}"/*.fastq "${RAW_DIR}"/*.fq)
if [[ ${#fastqs[@]} -eq 0 ]]; then
  echo "[ERROR] No FASTQ files found in ${RAW_DIR}"
  exit 1
fi

for fq in "${fastqs[@]}"; do
  base="$(basename "${fq}")"
  sample="${base}"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fq.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"
  sample="${sample%_1}"
  sample="${sample%_2}"
  out_prefix="${FILTERED_DIR}/${sample}"
  final_fq="${FILTERED_DIR}/${sample}.fastq.gz"

  if [[ -s "${final_fq}" ]]; then
    echo "[SKIP] ${final_fq} already exists."
    continue
  fi

  echo "[INFO] PRINSEQ filtering ${fq} -> ${final_fq}"
  tmp_fq="${FILTERED_DIR}/${sample}.tmp.fastq"
  if [[ "${fq}" == *.gz ]]; then
    gzip -dc "${fq}" > "${tmp_fq}"
  else
    cp "${fq}" "${tmp_fq}"
  fi

  ${PRINSEQ_CMD} \
    -fastq "${tmp_fq}" \
    -out_good "${out_prefix}" \
    -out_bad null \
    -min_len "${PRINSEQ_MIN_LEN}" \
    -min_qual_mean "${PRINSEQ_MIN_QUAL_MEAN}" \
    ${PRINSEQ_EXTRA_OPTS}

  rm -f "${tmp_fq}"
  if [[ ! -s "${out_prefix}.fastq" ]]; then
    echo "[ERROR] PRINSEQ did not create ${out_prefix}.fastq"
    exit 1
  fi
  gzip -f "${out_prefix}.fastq"
done

echo "[INFO] Building QIIME 2 single-end manifest..."
python3 - <<PY
from pathlib import Path
import csv
filtered = Path("${FILTERED_DIR}").resolve()
manifest = Path("${MANIFEST_TSV}")
fastqs = sorted(filtered.glob("*.fastq.gz"))
if not fastqs:
    raise SystemExit("No filtered fastq.gz files found in " + str(filtered))
with manifest.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample-id", "absolute-filepath", "direction"])
    for fq in fastqs:
        sample = fq.name.replace(".fastq.gz", "")
        writer.writerow([sample, str(fq), "forward"])
print(f"[INFO] Manifest written: {manifest}")
print(f"[INFO] Samples in manifest: {len(fastqs)}")
PY

echo "[INFO] Importing filtered reads into QIIME 2..."
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${MANIFEST_TSV}" \
  --output-path "${QIIME_DIR}/demux-prinseq.qza" \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "${QIIME_DIR}/demux-prinseq.qza" \
  --o-visualization "${QIIME_DIR}/demux-prinseq.qzv"

echo "[DONE] PRINSEQ-filtered reads imported into ${QIIME_DIR}/demux-prinseq.qza"
