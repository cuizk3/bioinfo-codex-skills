#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript deseq2_basic.R gene_counts.txt sample_info.csv outdir")
}

count_file <- args[1]
meta_file <- args[2]
outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading featureCounts output...")
raw <- read.delim(count_file, comment.char = "#", check.names = FALSE)
gene_ids <- raw$Geneid

count_cols <- grep(".sorted.bam$", colnames(raw), value = TRUE)
if (length(count_cols) == 0) {
  count_cols <- colnames(raw)[7:ncol(raw)]
}

counts <- raw[, count_cols, drop = FALSE]
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- sub(".sorted.bam$", "", colnames(counts))
rownames(counts) <- gene_ids

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
required <- c("sample", "group")
if (!all(required %in% colnames(meta))) {
  stop("metadata must contain columns: sample, group")
}

meta <- meta[match(colnames(counts), meta$sample), ]
if (any(is.na(meta$sample))) {
  stop("Some count matrix columns do not match metadata sample names.")
}
rownames(meta) <- meta$sample
meta$group <- factor(meta$group)

write.csv(counts, file.path(outdir, "count_matrix.csv"))
write.csv(meta, file.path(outdir, "metadata_matched.csv"))

message("Running DESeq2...")
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ group
)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(vst_mat, file.path(outdir, "vst_matrix.csv"))

# PCA
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14)
ggsave(file.path(outdir, "PCA.pdf"), p, width = 6, height = 5)
ggsave(file.path(outdir, "PCA.png"), p, width = 6, height = 5, dpi = 300)

# Sample correlation heatmap
pdf(file.path(outdir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(cor(vst_mat), annotation_col = meta[, "group", drop = FALSE])
dev.off()
png(file.path(outdir, "sample_correlation_heatmap.png"), width = 1800, height = 1500, res = 300)
pheatmap(cor(vst_mat), annotation_col = meta[, "group", drop = FALSE])
dev.off()

# Contrast: second group vs first group
groups <- levels(meta$group)
if (length(groups) < 2) {
  stop("At least two groups are required for differential expression.")
}
contrast <- c("group", groups[2], groups[1])
res <- results(dds, contrast = contrast)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

res_df <- res_df %>%
  mutate(
    status = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(res_df, file.path(outdir, paste0("DEG_", groups[2], "_vs_", groups[1], ".csv")), row.names = FALSE)

# Volcano plot
vol <- res_df %>%
  mutate(minus_log10_padj = -log10(padj),
         label = ifelse(rank(padj, ties.method = "first") <= 10, gene_id, NA))

p <- ggplot(vol, aes(log2FoldChange, minus_log10_padj, color = status)) +
  geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = label), max.overlaps = 20, na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    title = paste(groups[2], "vs", groups[1]),
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value"
  )
ggsave(file.path(outdir, "volcano.pdf"), p, width = 7, height = 6)
ggsave(file.path(outdir, "volcano.png"), p, width = 7, height = 6, dpi = 300)

# DEG heatmap
deg_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

if (length(deg_genes) >= 2) {
  mat <- vst_mat[deg_genes, , drop = FALSE]
  mat <- t(scale(t(mat)))
  pdf(file.path(outdir, "top50_DEG_heatmap.pdf"), width = 8, height = 10)
  pheatmap(mat, annotation_col = meta[, "group", drop = FALSE], show_rownames = TRUE)
  dev.off()
  png(file.path(outdir, "top50_DEG_heatmap.png"), width = 2100, height = 2700, res = 300)
  pheatmap(mat, annotation_col = meta[, "group", drop = FALSE], show_rownames = TRUE)
  dev.off()
}

# Summary
summary_file <- file.path(outdir, "analysis_summary.txt")
up_n <- sum(res_df$status == "Up", na.rm = TRUE)
down_n <- sum(res_df$status == "Down", na.rm = TRUE)
cat(
  "DESeq2 analysis summary\n",
  "Contrast: ", groups[2], " vs ", groups[1], "\n",
  "Up-regulated genes: ", up_n, "\n",
  "Down-regulated genes: ", down_n, "\n",
  "Note: Please combine DEG results with GO/KEGG/GSEA and known fish sex differentiation genes for biological interpretation.\n",
  file = summary_file,
  sep = ""
)

message("DESeq2 analysis completed: ", outdir)
