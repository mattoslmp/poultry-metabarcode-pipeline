#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${RAW_DIR}" "${METADATA_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/01_download_fastq_${BIOPROJECT}.log") 2>&1

ENA_REPORT="${METADATA_DIR}/ena_fastq_report.tsv"
ENA_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${BIOPROJECT}&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,experiment_alias,library_name,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&download=true"

echo "[INFO] BioProject: ${BIOPROJECT}"
echo "[INFO] Download method: ENA direct FASTQ links"
echo "[INFO] Downloading ENA FASTQ report..."

if ! curl -L --retry 5 --connect-timeout 30 "${ENA_URL}" -o "${ENA_REPORT}"; then
  if [[ -s "${ENA_REPORT}" ]]; then
    echo "[WARN] ENA query failed, using cached report: ${ENA_REPORT}"
  else
    echo "[ERROR] ENA query failed and no cached report exists: ${ENA_REPORT}"
    exit 1
  fi
fi

python3 - <<PY
import csv
from pathlib import Path

report = Path("${ENA_REPORT}")
runs_out = Path("${RUN_ACCESSIONS}")
meta_out = Path("${AUTO_METADATA_TSV}")
paths_out = Path("${METADATA_DIR}/ena_fastq_download_paths.tsv")

rows = []
with report.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        run = (row.get("run_accession") or "").strip()
        fastq = (row.get("fastq_ftp") or "").strip()
        if run and fastq:
            rows.append(row)

if not rows:
    raise SystemExit(f"No ENA FASTQ links found in {report}")

known_sites = {
    "SRR5975917": "Lung",
    "SRR5975918": "Trachea",
    "SRR5975919": "Lung",
    "SRR5975920": "Trachea",
    "SRR5975921": "Lung",
    "SRR5975922": "Trachea",
}

def infer_site(row):
    run = (row.get("run_accession") or "").strip()
    if run in known_sites:
        return known_sites[run]
    text = " ".join(str(v) for v in row.values()).lower()
    if any(x in text for x in ["trachea", "tracheal", "trach", "traqueia"]):
        return "Trachea"
    if any(x in text for x in ["lung", "lungs", "pulmao", "pulmão"]):
        return "Lung"
    return "Unknown"

runs_out.write_text("\n".join((r.get("run_accession") or "").strip() for r in rows) + "\n")

with meta_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample-id", "body_site", "sample_accession", "sample_alias", "sample_title"])
    for row in rows:
        writer.writerow([
            (row.get("run_accession") or "").strip(),
            infer_site(row),
            (row.get("sample_accession") or "").strip(),
            (row.get("sample_alias") or "").strip(),
            (row.get("sample_title") or row.get("experiment_alias") or row.get("library_name") or "").strip(),
        ])

with paths_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes"])
    for row in rows:
        writer.writerow([
            (row.get("run_accession") or "").strip(),
            (row.get("fastq_ftp") or "").strip(),
            (row.get("fastq_md5") or "").strip(),
            (row.get("fastq_bytes") or "").strip(),
        ])

print(f"[INFO] Runs found: {len(rows)}")
print(f"[INFO] Run list: {runs_out}")
print(f"[INFO] Auto metadata: {meta_out}")
print(f"[INFO] ENA download paths: {paths_out}")
PY

echo "[INFO] Downloading FASTQ files from ENA..."
tail -n +2 "${METADATA_DIR}/ena_fastq_download_paths.tsv" | while IFS=$'\t' read -r RUN FASTQFTP MD5 BYTES; do
  [[ -z "${RUN}" || -z "${FASTQFTP}" ]] && continue
  echo "============================================================"
  echo "[INFO] Run: ${RUN}"

  IFS=';' read -ra URLS <<< "${FASTQFTP}"
  for URL in "${URLS[@]}"; do
    [[ -z "${URL}" ]] && continue
    FILE="$(basename "${URL}")"
    OUT="${RAW_DIR}/${FILE}"

    if [[ -s "${OUT}" ]]; then
      echo "[SKIP] FASTQ already exists: ${OUT}"
      continue
    fi

    echo "[INFO] Downloading ${URL} -> ${OUT}"
    if ! wget -c --tries=10 --timeout=60 "https://${URL}" -O "${OUT}"; then
      wget -c --tries=10 --timeout=60 "ftp://${URL}" -O "${OUT}"
    fi
  done
done

if ! compgen -G "${RAW_DIR}/*.fastq.gz" >/dev/null; then
  echo "[ERROR] No FASTQ files were downloaded to ${RAW_DIR}"
  exit 1
fi

echo "[DONE] FASTQ files saved in ${RAW_DIR}"
ls -lh "${RAW_DIR}"/*.fastq.gz
