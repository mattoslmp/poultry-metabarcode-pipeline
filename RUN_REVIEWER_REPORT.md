# How to generate the reviewer sequencing report

After running the main FASTQ/QIIME 2 pipeline, generate the reviewer-ready sequencing report with:

```bash
conda activate qiime2-amplicon-2024.10
cd poultry-metabarcode-pipeline
python3 scripts/05_make_reviewer_report.py
```

The script writes:

```text
results/reviewer_report/reviewer_sequencing_summary.tsv
results/reviewer_report/reviewer_sequencing_summary.md
```

The report includes, per sample:

- raw FASTQ reads, if the raw FASTQ files are present in `data/raw_fastq`;
- PRINSEQ-filtered reads counted directly from `data/filtered_prinseq`;
- DADA2 denoising statistics exported from QIIME 2;
- non-rarefied feature-table reads and ASVs per sample;
- rarefied feature-table reads and ASVs per sample;
- retained percentage after rarefaction;
- whether the sample was included at the 7,503-read rarefaction depth.

Rarefaction is justified as a sequencing-depth normalization step before Shannon diversity estimation and Kruskal-Wallis testing.
