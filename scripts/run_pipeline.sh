#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

bash scripts/00_check_dependencies.sh
bash scripts/01_download_sra_prjna400142.sh
bash scripts/02_prinseq_filter_and_import.sh
bash scripts/03_qiime2_cutadapt_dada2_taxonomy.sh
python3 scripts/04_make_figures.py
python3 scripts/05_make_reviewer_report.py

echo "Pipeline complete. Results are in: ${PROJECT_DIR}/results"
