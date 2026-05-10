---
name: rnaseq
description: Run and review bulk RNA-seq workflows for fish transcriptomics, including QC, trimming, alignment, counting, DESeq2 differential expression, PCA, volcano plots, heatmaps, and biological interpretation.
---

# RNA-seq Skill

Use this skill when the user asks for:
- RNA-seq analysis
- nascent RNA-seq count analysis
- differential expression
- PCA, volcano plot, heatmap
- sex-biased gene expression in fish gonads
- dmrt1/foxl2/figla/cyp19a1a expression comparison

## Required inputs

- sample metadata CSV with columns:
  - `sample`
  - `group`
  - `fastq1`
  - `fastq2`
- reference genome FASTA
- gene annotation GTF
- or an existing gene count matrix and metadata

## Standard workflow

1. Validate metadata.
2. Run fastp for trimming and QC.
3. Run FastQC/MultiQC if available.
4. Build HISAT2 index if missing.
5. Align reads with HISAT2.
6. Sort/index BAM with samtools.
7. Count reads with featureCounts.
8. Run DESeq2.
9. Generate PCA, sample correlation heatmap, volcano plot and DEG heatmap.
10. Summarize top DEGs and sex-related genes.

## Biological interpretation rules

For fish sex differentiation:
- Male pathway genes: dmrt1, amh, amhr2, sox9, gsdf
- Female pathway genes: foxl2, figla, cyp19a1a, rspo1, wnt4
- Germ cell genes: vasa/ddx4, dazl, nanos3, piwil1

Do not claim causal mechanisms from expression alone. Say "associated with", "suggests", or "supports".
