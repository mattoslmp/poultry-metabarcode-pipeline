#!/usr/bin/env bash
# Main configuration for the PRJNA400142 poultry respiratory 16S pipeline.

# Project and compute settings
BIOPROJECT="${BIOPROJECT:-PRJNA400142}"
THREADS="${THREADS:-8}"

# PRINSEQ-lite filtering.
# The manuscript states reads >100 bp and mean Phred quality >= 30.
# Therefore the implemented minimum length is 101 and the mean quality cutoff is 30.
PRINSEQ_MIN_LEN="${PRINSEQ_MIN_LEN:-101}"
PRINSEQ_MIN_QUAL_MEAN="${PRINSEQ_MIN_QUAL_MEAN:-30}"
PRINSEQ_EXTRA_OPTS="${PRINSEQ_EXTRA_OPTS:-}"

# V4 16S primers 515F/806R and their reverse complements.
PRIMER_515F="${PRIMER_515F:-GTGYCAGCMGCCGCGGTAA}"
PRIMER_806R="${PRIMER_806R:-GGACTACNVGGGTWTCTAAT}"
PRIMER_515F_RC="${PRIMER_515F_RC:-TTACCGCGGCKGCTGRCAC}"
PRIMER_806R_RC="${PRIMER_806R_RC:-ATTAGAWACCCBNGTAGTCC}"
CUTADAPT_OVERLAP="${CUTADAPT_OVERLAP:-15}"
CUTADAPT_ERROR_RATE="${CUTADAPT_ERROR_RATE:-0.1}"

# DADA2 single-end settings.
DADA2_TRIM_LEFT="${DADA2_TRIM_LEFT:-0}"
DADA2_TRUNC_LEN="${DADA2_TRUNC_LEN:-0}"
DADA2_CHIMERA_METHOD="${DADA2_CHIMERA_METHOD:-consensus}"

# Successful rarefaction depth used for the reviewer response.
# This depth retains samples from both anatomical sites after filtering/denoising.
RAREFY_DEPTH="${RAREFY_DEPTH:-3743}"

# Classifier path. Place the classifier here before running taxonomy.
PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CLASSIFIER_QZA="${CLASSIFIER_QZA:-${PROJECT_DIR}/reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza}"
CLASSIFIER_URL="${CLASSIFIER_URL:-}"

# Directory layout.
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/data/raw_fastq}"
FILTERED_DIR="${FILTERED_DIR:-${PROJECT_DIR}/data/filtered_prinseq}"
METADATA_DIR="${METADATA_DIR:-${PROJECT_DIR}/metadata}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
QIIME_DIR="${QIIME_DIR:-${RESULTS_DIR}/qiime}"
EXPORTED_DIR="${EXPORTED_DIR:-${RESULTS_DIR}/exported}"
FIGURES_DIR="${FIGURES_DIR:-${RESULTS_DIR}/figures}"
LOG_DIR="${LOG_DIR:-${RESULTS_DIR}/logs}"

RUNINFO_CSV="${RUNINFO_CSV:-${METADATA_DIR}/${BIOPROJECT}_RunInfo.csv}"
RUN_ACCESSIONS="${RUN_ACCESSIONS:-${METADATA_DIR}/run_accessions.txt}"
AUTO_METADATA_TSV="${AUTO_METADATA_TSV:-${METADATA_DIR}/sample-metadata.auto.tsv}"
MANUAL_METADATA_TSV="${MANUAL_METADATA_TSV:-${METADATA_DIR}/sample-metadata.tsv}"
MANIFEST_TSV="${MANIFEST_TSV:-${METADATA_DIR}/manifest-single-end.tsv}"
