#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${PROJECT_DIR}/config/pipeline_config.sh"

mkdir -p "${RAW_DIR}" "${FILTERED_DIR}" "${METADATA_DIR}" "${QIIME_DIR}" "${EXPORTED_DIR}" "${FIGURES_DIR}" "${LOG_DIR}" "${PROJECT_DIR}/reference"

need_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "[MISSING] ${cmd}"
    return 1
  fi
  echo "[OK] ${cmd}: $(command -v "${cmd}")"
}

missing=0
for cmd in qiime python3 prefetch fasterq-dump esearch efetch gzip awk sed biom; do
  need_cmd "${cmd}" || missing=1
done

if command -v prinseq-lite.pl >/dev/null 2>&1; then
  echo "[OK] prinseq-lite.pl: $(command -v prinseq-lite.pl)"
elif command -v prinseq-lite >/dev/null 2>&1; then
  echo "[OK] prinseq-lite: $(command -v prinseq-lite)"\else
  echo "[MISSING] prinseq-lite.pl or prinseq-lite"
  missing=1
fi

python3 - <<'PY'
import importlib
mods = ["pandas", "numpy", "scipy", "matplotlib"]
missing = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception:
        missing.append(m)
if missing:
    raise SystemExit("Missing Python modules: " + ", ".join(missing))
print("[OK] Python modules: pandas, numpy, scipy, matplotlib")
try:
    import matplotlib_venn
    print("[OK] Python module: matplotlib_venn")
except Exception:
    print("[WARN] matplotlib_venn not found; the figure script will use a fallback Venn drawing.")
PY

if [[ ! -s "${CLASSIFIER_QZA}" ]]; then
  if [[ -n "${CLASSIFIER_URL}" ]]; then
    echo "[INFO] Classifier not found. Downloading from CLASSIFIER_URL..."
    mkdir -p "$(dirname "${CLASSIFIER_QZA}")"
    wget -O "${CLASSIFIER_QZA}" "${CLASSIFIER_URL}"
  else
    echo "[MISSING] Classifier not found: ${CLASSIFIER_QZA}"
    echo "Place the Greengenes2 2024.09 V4 classifier there or set CLASSIFIER_QZA=/path/to/classifier.qza"
    missing=1
  fi
else
  echo "[OK] Classifier: ${CLASSIFIER_QZA}"
fi

if [[ "${missing}" -ne 0 ]]; then
  echo
  echo "Some dependencies are missing. Suggested install command inside the QIIME 2 2024.10 environment:"
  echo "conda install -c conda-forge -c bioconda sra-tools entrez-direct prinseq-lite biom-format pandas numpy scipy matplotlib matplotlib-venn -y"
  exit 1
fi

echo "All required dependencies were found."
