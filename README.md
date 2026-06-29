# Poultry respiratory 16S metabarcode pipeline — PRJNA400142

Reproducible pipeline to download the raw FASTQ data from NCBI SRA BioProject **PRJNA400142**, filter reads with **PRINSEQ-lite**, process V4 16S single-end reads in **QIIME 2 2024.10**, classify ASVs with a **Greengenes2 2024.09 V4 Naive Bayes classifier**, and generate the article-style figures for lung vs. trachea respiratory microbiomes.

The pipeline was written to reproduce the methodology described in the manuscript *Bacterial community of the respiratory tract of clinically healthy broilers*.

## Main outputs

After a successful run, the main outputs will be written to `results/`:

- `results/qiime/`: QIIME 2 artifacts and visualizations.
- `results/exported/`: exported feature tables, taxonomy, alpha diversity vectors, and stats.
- `results/figures/Figure_3_panel.png`, `.pdf`, `.svg`: article-style panel with:
  - A: stacked bar chart at class level, QIIME taxonomy level 3;
  - B: stacked bar chart at genus level, QIIME taxonomy level 6;
  - C: Venn diagram of genera detected in lung and trachea;
  - D: genus-level pseudo-log10 heatmap by location.
- `results/figures/shannon_alpha_diversity.png`, `.pdf`, `.svg`: rarefied Shannon diversity comparison.
- `results/exported/shannon_kruskal_wallis.tsv`: Kruskal–Wallis H statistic and p-value.

## Repository structure

```text
config/
  pipeline_config.sh              # Main editable parameters
metadata/
  sample-metadata.tsv             # Manual metadata template for QIIME 2
scripts/
  00_check_dependencies.sh
  01_download_sra_prjna400142.sh
  02_prinseq_filter_and_import.sh
  03_qiime2_cutadapt_dada2_taxonomy.sh
  04_make_figures.py
  run_pipeline.sh
METHOD.md                         # Methodology text saved separately
README.md
```

## 1. Installation

### Option A — use an existing QIIME 2 environment

If you already have QIIME 2 2024.10 installed, activate it before running the pipeline:

```bash
conda activate qiime2-amplicon-2024.10
```

Then install the additional command-line and Python utilities if they are missing:

```bash
conda install -c conda-forge -c bioconda sra-tools entrez-direct prinseq-lite biom-format pandas numpy scipy matplotlib matplotlib-venn -y
```

### Option B — create QIIME 2 2024.10 from the official environment file

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py39-linux-conda.yml
conda env create -n qiime2-amplicon-2024.10 --file qiime2-amplicon-2024.10-py39-linux-conda.yml
conda activate qiime2-amplicon-2024.10
conda install -c conda-forge -c bioconda sra-tools entrez-direct prinseq-lite biom-format pandas numpy scipy matplotlib matplotlib-venn -y
```

## 2. Clone and enter the repository

```bash
git clone https://github.com/mattoslmp/poultry-metabarcode-pipeline.git
cd poultry-metabarcode-pipeline
chmod +x scripts/*.sh
```

## 3. Prepare the Greengenes2 classifier

The taxonomy step expects a pretrained V4 classifier compatible with **Greengenes2 2024.09** and **scikit-learn 1.4.2**.

Place it here:

```text
reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza
```

or edit `config/pipeline_config.sh` and set:

```bash
CLASSIFIER_QZA=/absolute/path/to/your/classifier.qza
```

If you have a stable download URL for the classifier, set:

```bash
export CLASSIFIER_URL="https://.../classifier.qza"
```

and the pipeline will download it automatically.

## 4. Run the complete pipeline

```bash
conda activate qiime2-amplicon-2024.10
bash scripts/run_pipeline.sh
```

## 5. Run step by step

```bash
bash scripts/00_check_dependencies.sh
bash scripts/01_download_sra_prjna400142.sh
bash scripts/02_prinseq_filter_and_import.sh
bash scripts/03_qiime2_cutadapt_dada2_taxonomy.sh
python scripts/04_make_figures.py
```

## 6. Metadata checking

The download script creates:

```text
metadata/PRJNA400142_RunInfo.csv
metadata/sample-metadata.auto.tsv
metadata/run_accessions.txt
```

The pipeline uses `metadata/sample-metadata.tsv` if it exists and contains real sample rows. Otherwise, it falls back to `metadata/sample-metadata.auto.tsv`.

Before publication-quality analysis, manually check that the `body_site` column correctly labels each SRA run as `Lung` or `Trachea`.

## 7. Methodological target implemented

This repository implements the following workflow:

1. Download raw single-end reads from NCBI SRA BioProject PRJNA400142.
2. Filter raw reads with PRINSEQ-lite, retaining reads longer than 100 bp and with mean Phred quality score >= 30.
3. Import filtered single-end reads into QIIME 2.
4. Remove 515F/806R V4 primers and reverse complements using q2-cutadapt with 5-prime anchored matching and minimum 15-nt overlap.
5. Denoise with DADA2 single-end mode using `trim-left = 0`, `trunc-len = 0`, and consensus chimera removal.
6. Classify representative ASVs using `classify-sklearn` and the Greengenes2 2024.09 V4 classifier.
7. Collapse feature tables to class (L3) and genus (L6) ranks.
8. Generate stacked bar charts for lung vs. trachea at class and genus levels.
9. Generate genus-level Venn and heatmap figures.
10. Rarefy to 7,503 reads per sample and calculate Shannon diversity.
11. Compare lung vs. trachea Shannon diversity with a two-sided Kruskal–Wallis test.

## 8. Notes

- This is a low-biomass respiratory microbiome dataset. Interpret environmental taxa cautiously.
- The pipeline does not discard untrimmed reads during cutadapt unless you edit the config and add that option manually.
- The default minimum length is `101` because the written method says reads `>100 bp`.
- The PRINSEQ filter uses `-min_qual_mean 30`; if you need stricter per-base filtering, add extra PRINSEQ options in `config/pipeline_config.sh`.

## Citation and data availability

Raw sequence data: NCBI SRA BioProject **PRJNA400142**.

Pipeline repository: https://github.com/mattoslmp/poultry-metabarcode-pipeline
