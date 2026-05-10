#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript enrichment_clusterprofiler.R DEG.csv gene2go.tsv gene2kegg.tsv outdir")
}

deg_file <- args[1]
gene2go_file <- args[2]
gene2kegg_file <- args[3]
outdir <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
required <- c("gene_id", "log2FoldChange", "padj")
if (!all(required %in% colnames(deg))) {
  stop("DEG table must contain gene_id, log2FoldChange and padj.")
}

sig <- deg %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

if (nrow(sig) < 5) {
  warning("Fewer than 5 significant DEGs. Enrichment may be unstable.")
}

# GO enrichment using custom TERM2GENE
if (file.exists(gene2go_file)) {
  gene2go <- read.delim(gene2go_file, header = TRUE, stringsAsFactors = FALSE)
  # required columns: go_id, gene_id
  if (!all(c("go_id", "gene_id") %in% colnames(gene2go))) {
    stop("gene2go.tsv must contain columns: go_id, gene_id")
  }

  ego <- enricher(
    gene = sig$gene_id,
    universe = deg$gene_id,
    TERM2GENE = gene2go[, c("go_id", "gene_id")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  go_df <- as.data.frame(ego)
  write.csv(go_df, file.path(outdir, "GO_enrichment.csv"), row.names = FALSE)

  if (nrow(go_df) > 0) {
    p <- dotplot(ego, showCategory = 20) + ggtitle("GO enrichment")
    ggsave(file.path(outdir, "GO_dotplot.pdf"), p, width = 8, height = 6)
    ggsave(file.path(outdir, "GO_dotplot.png"), p, width = 8, height = 6, dpi = 300)
  }
}

# KEGG-like enrichment using custom TERM2GENE
if (file.exists(gene2kegg_file)) {
  gene2kegg <- read.delim(gene2kegg_file, header = TRUE, stringsAsFactors = FALSE)
  # required columns: pathway_id, gene_id
  if (!all(c("pathway_id", "gene_id") %in% colnames(gene2kegg))) {
    stop("gene2kegg.tsv must contain columns: pathway_id, gene_id")
  }

  ekegg <- enricher(
    gene = sig$gene_id,
    universe = deg$gene_id,
    TERM2GENE = gene2kegg[, c("pathway_id", "gene_id")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  kegg_df <- as.data.frame(ekegg)
  write.csv(kegg_df, file.path(outdir, "KEGG_enrichment.csv"), row.names = FALSE)

  if (nrow(kegg_df) > 0) {
    p <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG enrichment")
    ggsave(file.path(outdir, "KEGG_dotplot.pdf"), p, width = 8, height = 6)
    ggsave(file.path(outdir, "KEGG_dotplot.png"), p, width = 8, height = 6, dpi = 300)
  }
}

# GSEA using ranked log2FC and custom GO
if (file.exists(gene2go_file)) {
  gene2go <- read.delim(gene2go_file, header = TRUE, stringsAsFactors = FALSE)
  gene_list <- deg$log2FoldChange
  names(gene_list) <- deg$gene_id
  gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)

  gsea_go <- GSEA(
    geneList = gene_list,
    TERM2GENE = gene2go[, c("go_id", "gene_id")],
    pvalueCutoff = 0.05,
    verbose = FALSE
  )

  gsea_df <- as.data.frame(gsea_go)
  write.csv(gsea_df, file.path(outdir, "GSEA_GO.csv"), row.names = FALSE)

  if (nrow(gsea_df) > 0) {
    p <- dotplot(gsea_go, showCategory = 20) + ggtitle("GSEA GO")
    ggsave(file.path(outdir, "GSEA_GO_dotplot.pdf"), p, width = 8, height = 6)
    ggsave(file.path(outdir, "GSEA_GO_dotplot.png"), p, width = 8, height = 6, dpi = 300)
  }
}

cat("Enrichment analysis completed: ", outdir, "\n")
