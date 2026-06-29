#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${RAW_DIR}" "${METADATA_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_DIR}/01_download_sra_${BIOPROJECT}.log") 2>&1

echo "[INFO] BioProject: ${BIOPROJECT}"
echo "[INFO] Querying NCBI SRA run info..."

esearch -db sra -query "${BIOPROJECT}" | efetch -format runinfo > "${RUNINFO_CSV}"

if [[ ! -s "${RUNINFO_CSV}" ]]; then
  echo "[ERROR] Empty run-info table: ${RUNINFO_CSV}"
  exit 1
fi

python3 - <<PY
import csv
from pathlib import Path
runinfo = Path("${RUNINFO_CSV}")
runs_out = Path("${RUN_ACCESSIONS}")
meta_out = Path("${AUTO_METADATA_TSV}")
rows = []
with runinfo.open(newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        run = (row.get("Run") or "").strip()
        if run:
            rows.append(row)
if not rows:
    raise SystemExit("No SRA runs found in " + str(runinfo))
runs_out.write_text("\n".join((r.get("Run") or "").strip() for r in rows if (r.get("Run") or "").strip()) + "\n")

def infer_site(row):
    text = " ".join(str(row.get(k, "")) for k in [
        "Run", "SampleName", "LibraryName", "LibraryLayout", "Title", "source_name", "tissue", "BioSample", "Experiment"
    ]).lower()
    if any(x in text for x in ["trachea", "tracheal", "traqueia"]):
        return "Trachea"
    if any(x in text for x in ["lung", "lungs", "pulmao", "pulmão"]):
        return "Lung"
    return "Unknown"

with meta_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample-id", "body_site", "biosample", "sra_title"])
    for row in rows:
        writer.writerow([
            (row.get("Run") or "").strip(),
            infer_site(row),
            (row.get("BioSample") or "").strip(),
            (row.get("Title") or row.get("LibraryName") or "").strip(),
        ])
print(f"[INFO] Runs found: {len(rows)}")
print(f"[INFO] Run list: {runs_out}")
print(f"[INFO] Auto metadata: {meta_out}")
PY

echo "[INFO] Downloading runs with prefetch and fasterq-dump..."
while read -r run; do
  [[ -z "${run}" ]] && continue
  echo "[INFO] Processing ${run}"
  if compgen -G "${RAW_DIR}/${run}*.fastq.gz" >/dev/null; then
    echo "[SKIP] FASTQ already exists for ${run}"
    continue
  fi
  prefetch --max-size "${MAX_SRA_SIZE}" "${run}"
  fasterq-dump --threads "${THREADS}" --outdir "${RAW_DIR}" --split-3 "${run}"
  gzip -f "${RAW_DIR}/${run}"*.fastq
done < "${RUN_ACCESSIONS}"

echo "[DONE] FASTQ files saved in ${RAW_DIR}"
