#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_gene_expression.R vst_matrix.csv metadata.csv gene_id outdir")
}

expr_file <- args[1]
meta_file <- args[2]
gene_id <- args[3]
outdir <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

if (!gene_id %in% rownames(expr)) {
  stop("gene_id not found in expression matrix: ", gene_id)
}

df <- data.frame(
  sample = colnames(expr),
  expression = as.numeric(expr[gene_id, ])
) %>%
  left_join(meta, by = "sample")

p <- ggplot(df, aes(x = group, y = expression)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(title = gene_id, x = "Group", y = "VST expression")

ggsave(file.path(outdir, paste0(gene_id, "_expression.pdf")), p, width = 5, height = 4)
ggsave(file.path(outdir, paste0(gene_id, "_expression.png")), p, width = 5, height = 4, dpi = 300)
