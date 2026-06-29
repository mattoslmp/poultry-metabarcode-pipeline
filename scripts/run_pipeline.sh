#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

export TMPDIR="${TMPDIR:-${PROJECT_DIR}/tmp_qiime}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

# Successful workflow used for the reviewer response.
bash scripts/00_check_dependencies.sh
bash scripts/01_download_sra_prjna400142.sh
bash scripts/02_prinseq_filter_and_import.sh
bash scripts/03a_cutadapt_dada2_only.sh
bash scripts/03b_resume_after_dada2_lowmem_taxonomy.sh
bash scripts/03c_resume_after_taxonomy_no_ordination.sh
python3 scripts/05_make_reviewer_report.py
python3 scripts/04_make_figures.py

mkdir -p results/figures/exact_from_article_script
python3 scripts/make_all_panels_exact.py \
  --table results/qiime/table.qza \
  --taxonomy results/qiime/taxonomy.qza \
  --metadata results/exported/sample-metadata.used.tsv \
  --group-col body_site \
  --group-a Lung \
  --group-b Trachea \
  --outdir results/figures/exact_from_article_script

echo "Pipeline complete. Results are in: ${PROJECT_DIR}/results"
