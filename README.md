# Poultry respiratory 16S metabarcode pipeline — PRJNA400142

This repository contains the reproducible workflow used to document the 16S rRNA amplicon sequencing analysis for BioProject **PRJNA400142**. It was organized for the reviewer response and includes only the steps that were successfully executed during the final reanalysis.

The workflow downloads the single-end FASTQ files directly from **ENA FASTQ links**, filters reads with **PRINSEQ-lite**, processes the V4 16S reads with **QIIME 2 2024.10** and **DADA2**, classifies ASVs with a **Greengenes2 2024.09 V4 Naive Bayes classifier**, exports reviewer-ready read/ASV tables, and generates Figure 3-style taxonomic panels.

## Main reviewer outputs

After the successful workflow finishes, the files most relevant to the reviewer response are:

```text
results/reviewer_report/reviewer_sequencing_summary.tsv
results/reviewer_report/reviewer_sequencing_summary.md
results/exported/feature_table.tsv
results/exported/rarefied_table.tsv
results/exported/denoising_stats/stats.tsv
results/exported/taxonomy/taxonomy.tsv
results/exported/taxa_level3.tsv
results/exported/taxa_level6.tsv
results/exported/shannon/alpha-diversity.tsv
results/figures/Figure_3_panel.png
results/figures/Figure_3_panel.pdf
results/figures/Figure_3_panel.svg
results/figures/exact_from_article_script/fig_ABCD_A4.png
results/figures/exact_from_article_script/fig_ABCD_A4.svg
```

The key reviewer table is:

```text
results/reviewer_report/reviewer_sequencing_summary.tsv
```

It reports, per sample, raw FASTQ reads, PRINSEQ-filtered reads, DADA2-filtered reads, denoised reads, non-chimeric reads, ASVs per sample, rarefied reads, ASVs retained after rarefaction, the rarefaction depth, retained percentage, and whether the sample was included in the rarefied alpha-diversity analysis.

## Repository structure

```text
config/
  pipeline_config.sh
metadata/
  sample-metadata.tsv
scripts/
  00_check_dependencies.sh
  01_download_sra_prjna400142.sh
  02_prinseq_filter_and_import.sh
  03a_cutadapt_dada2_only.sh
  03b_resume_after_dada2_lowmem_taxonomy.sh
  03c_resume_after_taxonomy_no_ordination.sh
  04_make_figures.py
  05_make_reviewer_report.py
  06_patch_exact_panel_script.sh
  make_all_panels_exact.py
  run_pipeline.sh
METHOD.md
REVIEWER_RESPONSE.md
RUN_REVIEWER_REPORT.md
README.md
```

## Successful workflow used in the final reanalysis

### 1. Activate QIIME 2 and set paths

```bash
cd ~/Benito/poultry-metabarcode-pipeline
conda activate qiime2-amplicon-2024.10
export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r
```

### 2. Install the required utilities

```bash
conda install -n qiime2-amplicon-2024.10 \
  -c bioconda -c conda-forge \
  prinseq matplotlib-venn biom-format pandas numpy scipy matplotlib wget curl -y
```

Check PRINSEQ:

```bash
which prinseq-lite.pl || find "$CONDA_PREFIX" -iname '*prinseq*' 2>/dev/null | head
```

### 3. Download the FASTQ files from ENA

The successful download path used ENA direct FASTQ links, not `prefetch` or `fasterq-dump`.

```bash
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

### 4. Use the checked sample metadata

```bash
cat > metadata/sample-metadata.tsv <<'EOF'
sample-id	body_site
SRR5975917	Lung
SRR5975918	Trachea
SRR5975919	Lung
SRR5975920	Trachea
SRR5975921	Lung
SRR5975922	Trachea
EOF
```

### 5. Download the Greengenes2 V4 classifier

```bash
mkdir -p reference
wget -c \
  https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.v4.nb.qza \
  -O reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza

qiime tools peek reference/gg2-2024.09-v4-classifier-sklearn-1.4.2.qza
```

The expected artifact type is:

```text
FeatureData[TaxonomicClassifier]
```

### 6. Create a QIIME temporary directory outside `/tmp`

This avoids failures caused by a full system `/tmp` directory.

```bash
mkdir -p tmp_qiime
rm -rf tmp_qiime/*
export TMPDIR="$PWD/tmp_qiime"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
```

### 7. Run PRINSEQ filtering and QIIME import

```bash
bash scripts/02_prinseq_filter_and_import.sh
```

This step applies the manuscript filter:

```text
PRINSEQ-lite: reads >100 bp; mean Phred quality >=30
```

### 8. Run primer trimming and DADA2 only

```bash
bash scripts/03a_cutadapt_dada2_only.sh
```

This generates:

```text
results/qiime/table.qza
results/qiime/rep-seqs.qza
results/qiime/denoising-stats.qza
```

### 9. Run low-memory taxonomy classification

```bash
export TAXONOMY_JOBS=1
export READS_PER_BATCH=50
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

bash scripts/03b_resume_after_dada2_lowmem_taxonomy.sh
```

This generates:

```text
results/qiime/taxonomy.qza
results/qiime/taxa_level3.qza
results/qiime/taxa_level6.qza
results/qiime/taxa-bar-plots.qzv
```

### 10. Rarefy and export results without beta ordination

The successful reviewer workflow used direct rarefaction and Shannon calculation, avoiding beta-diversity ordination because it was not required for the reviewer request.

```bash
export RAREFY_DEPTH=3743
bash scripts/03c_resume_after_taxonomy_no_ordination.sh
```

The rarefaction depth of **3,743 reads per sample** was selected because it retained samples from both anatomical sites after filtering, denoising and chimera removal.

### 11. Generate the reviewer table

```bash
python3 scripts/05_make_reviewer_report.py
```

Main output:

```text
results/reviewer_report/reviewer_sequencing_summary.tsv
```

### 12. Generate figures

```bash
python3 scripts/04_make_figures.py
```

To generate the Figure 3 layout using the plotting script supplied with the manuscript data package:

```bash
bash scripts/06_patch_exact_panel_script.sh
mkdir -p results/figures/exact_from_article_script

python3 scripts/make_all_panels_exact.py \
  --table results/qiime/table.qza \
  --taxonomy results/qiime/taxonomy.qza \
  --metadata results/exported/sample-metadata.used.tsv \
  --group-col body_site \
  --group-a Lung \
  --group-b Trachea \
  --outdir results/figures/exact_from_article_script
```

## One-command workflow

After installing dependencies and downloading the Greengenes2 classifier, the successful workflow can be run with:

```bash
bash scripts/run_pipeline.sh
```

## Notes for the reviewer response

The manuscript reports **358,104 fragments** as the initial sequencing output. The table reports the reads retained after the documented filtering, denoising and chimera-removal steps. These are different stages of the same sequencing workflow and should not be interpreted as the same count.


## Associated manuscript

**Title:** Bacterial community of the respiratory tract of clinically healthy broilers

**Authors:** Soares BD; Kunert-Filho HC; Grassotti TT; Gazal LES; Carvalho D; Pereira LM; Oliveira RR; Borges KA; Furian TQ; Nickel VS; Brito KCT; Kobayashi RKT; Destri GAD; Otutumi LK; Brito BG.

**Journal:** Brazilian Journal of Poultry Science / Revista Brasileira de Ciência Avícola.

**Status:** Manuscript submitted and currently under evaluation.

## Data availability

Raw sequence data: NCBI/ENA BioProject **PRJNA400142**.

Pipeline repository: https://github.com/mattoslmp/poultry-metabarcode-pipeline
