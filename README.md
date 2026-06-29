# Poultry respiratory 16S metabarcode pipeline — PRJNA400142

Reproducible pipeline for BioProject **PRJNA400142**. The pipeline downloads the raw single-end FASTQ files directly from **ENA FASTQ links** to avoid SSL/certificate problems observed with older `prefetch/fasterq-dump` installations, filters reads with **PRINSEQ-lite**, processes V4 16S reads with **QIIME 2 2024.10**, classifies ASVs with a **Greengenes2 2024.09 V4 Naive Bayes classifier**, generates the manuscript-style Figure 3 panels, and writes reviewer-ready sequencing summary tables.

The pipeline was written to reproduce and document the sequencing analysis described in the manuscript *Bacterial community of the respiratory tract of clinically healthy broilers*.

## Main outputs needed for reviewers

After a successful run, the key reviewer files are:

```text
results/reviewer_report/reviewer_sequencing_summary.tsv
results/reviewer_report/reviewer_sequencing_summary.md
results/exported/shannon_kruskal_wallis.tsv
results/exported/feature_table.tsv
results/exported/rarefied_table.tsv
results/exported/denoising_stats/stats.tsv
results/figures/Figure_3_panel.png
results/figures/Figure_3_panel.pdf
results/figures/Figure_3_panel.svg
results/figures/shannon_alpha_diversity.png
results/figures/shannon_alpha_diversity.pdf
results/figures/shannon_alpha_diversity.svg
```

The file `results/reviewer_report/reviewer_sequencing_summary.tsv` reports, per sample, the raw FASTQ counts, PRINSEQ-filtered reads, DADA2/QIIME 2 counts, non-rarefied reads and ASVs, rarefied reads and ASVs, retained percentage after rarefaction, and whether the sample was included at the rarefaction depth of 7,503 sequences.

## Repository structure

```text
config/
  pipeline_config.sh
metadata/
  sample-metadata.tsv
scripts/
  00_check_dependencies.sh
  01_download_sra_prjna400142.sh      # ENA FASTQ download for PRJNA400142
  02_prinseq_filter_and_import.sh
  03_qiime2_cutadapt_dada2_taxonomy.sh
  04_make_figures.py
  05_make_reviewer_report.py
  run_pipeline.sh
METHOD.md
REVIEWER_RESPONSE.md
RUN_REVIEWER_REPORT.md
README.md
```

## 1. Installation

Activate your QIIME 2 2024.10 environment:

```bash
conda activate qiime2-amplicon-2024.10
export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r
```

Install the extra tools required by this pipeline:

```bash
conda install -n qiime2-amplicon-2024.10 \
  -c bioconda -c conda-forge \
  prinseq matplotlib-venn biom-format pandas numpy scipy matplotlib wget curl -y
```

Check that PRINSEQ is available:

```bash
which prinseq-lite.pl || find "$CONDA_PREFIX" -iname '*prinseq*' 2>/dev/null | head
```

## 2. Clone or update the repository

```bash
git clone https://github.com/mattoslmp/poultry-metabarcode-pipeline.git
cd poultry-metabarcode-pipeline
chmod +x scripts/*.sh
```

If the repository already exists locally:

```bash
cd ~/Benito/poultry-metabarcode-pipeline
git pull origin main
chmod +x scripts/*.sh
```

If `git pull` is blocked by local changes, save them first:

```bash
git status --short
mkdir -p backup_local_changes
cp scripts/03_qiime2_cutadapt_dada2_taxonomy.sh backup_local_changes/03_qiime2_local_backup.sh 2>/dev/null || true
git stash push -m "backup local changes before update"
git pull origin main
chmod +x scripts/*.sh
```

## 3. Download FASTQ files from ENA

The download step now uses ENA direct FASTQ links and does not require `prefetch` or `fasterq-dump`:

```bash
cd ~/Benito/poultry-metabarcode-pipeline
conda activate qiime2-amplicon-2024.10
export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r

bash scripts/01_download_sra_prjna400142.sh
```

Expected FASTQ files:

```text
data/raw_fastq/SRR5975917.fastq.gz
data/raw_fastq/SRR5975918.fastq.gz
data/raw_fastq/SRR5975919.fastq.gz
data/raw_fastq/SRR5975920.fastq.gz
data/raw_fastq/SRR5975921.fastq.gz
data/raw_fastq/SRR5975922.fastq.gz
```

Check the download:

```bash
ls -lh data/raw_fastq/
cat metadata/run_accessions.txt
cat metadata/sample-metadata.auto.tsv
```

## 4. Prepare the Greengenes2 classifier

The taxonomy step expects the Greengenes2 2024.09 V4 classifier here:

```text
reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza
```

Download it:

```bash
mkdir -p reference
wget -c \
  https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.v4.nb.qza \
  -O reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza

qiime tools peek reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza
```

The expected QIIME 2 type is `FeatureData[TaxonomicClassifier]`.

## 5. Run the full analysis

Use `&&` so the pipeline stops if one step fails:

```bash
bash scripts/02_prinseq_filter_and_import.sh && \
bash scripts/03_qiime2_cutadapt_dada2_taxonomy.sh && \
python3 scripts/04_make_figures.py && \
python3 scripts/05_make_reviewer_report.py
```

Alternatively, once all dependencies and the classifier are present:

```bash
bash scripts/run_pipeline.sh
```

## 6. Files to send or use for the reviewer response

Use these files for the reviewer response and revised manuscript:

```bash
ls -lh results/reviewer_report/
cat results/reviewer_report/reviewer_sequencing_summary.tsv
cat results/reviewer_report/reviewer_sequencing_summary.md
ls -lh results/figures/
```

Important files:

- `results/reviewer_report/reviewer_sequencing_summary.tsv`: table requested by reviewers.
- `results/reviewer_report/reviewer_sequencing_summary.md`: reviewer-friendly markdown report.
- `results/exported/shannon_kruskal_wallis.tsv`: H statistic and p-value.
- `results/figures/Figure_3_panel.*`: manuscript-style Figure 3.
- `results/figures/figure3A_class_barplot.*`: Figure 3A.
- `results/figures/figure3B_genus_barplot.*`: Figure 3B.
- `results/figures/figure3C_venn.*`: Figure 3C.
- `results/figures/figure3D_genus_heatmap.*`: Figure 3D.
- `results/figures/shannon_alpha_diversity.*`: rarefied Shannon diversity figure.

## 7. Methodological target implemented

This repository implements the following workflow:

1. Download raw single-end FASTQ reads for BioProject PRJNA400142 using ENA direct FASTQ links.
2. Filter raw reads with PRINSEQ-lite, retaining reads longer than 100 bp and with mean Phred quality score >= 30.
3. Import filtered single-end reads into QIIME 2.
4. Remove 515F/806R V4 primers and reverse complements using q2-cutadapt with 5-prime anchored matching and minimum 15-nt overlap.
5. Denoise with DADA2 single-end mode using `trim-left = 0`, `trunc-len = 0`, and consensus chimera removal.
6. Classify representative ASVs using `classify-sklearn` and the Greengenes2 2024.09 V4 classifier.
7. Collapse feature tables to class (L3) and genus (L6) ranks.
8. Generate stacked bar charts for lung vs. trachea at class and genus levels.
9. Generate genus-level Venn and heatmap figures.
10. Rarefy to 7,503 reads per sample and calculate Shannon diversity.
11. Compare lung vs. trachea Shannon diversity with a two-sided Kruskal-Wallis test.
12. Generate reviewer-ready read-count and ASV-count summary tables.

## 8. Rarefaction justification

Rarefaction to 7,503 sequences per sample is used to normalize sequencing depth across samples before diversity estimation. This standardization reduces sampling-depth bias and allows Shannon diversity and Kruskal-Wallis comparisons to be performed at a common sequencing depth across lung and trachea samples. Samples with fewer reads than the selected depth are excluded by the QIIME 2 core-metrics workflow.

## Citation and data availability

Raw sequence data: BioProject **PRJNA400142**.

Pipeline repository: https://github.com/mattoslmp/poultry-metabarcode-pipeline
