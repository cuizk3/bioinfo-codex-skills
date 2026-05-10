#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript wgcna_basic.R vst_matrix.csv traits.csv outdir")
}

expr_file <- args[1]
trait_file <- args[2]
outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
traits <- read.csv(trait_file, row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(expr), rownames(traits))
if (length(common_samples) < 6) {
  stop("WGCNA usually needs more samples. Fewer than 6 matched samples found.")
}
expr <- expr[, common_samples, drop = FALSE]
traits <- traits[common_samples, , drop = FALSE]

# genes as columns for WGCNA
datExpr <- as.data.frame(t(expr))

# filter low variance genes
vars <- apply(datExpr, 2, var)
keep <- vars >= quantile(vars, 0.75, na.rm = TRUE)
datExpr <- datExpr[, keep]

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(outdir, "sample_clustering.pdf"), width = 8, height = 5)
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "")
dev.off()

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file.path(outdir, "soft_threshold.pdf"), width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold power", ylab = "Scale Free Topology Model Fit",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.8)
abline(h = 0.8, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold power", ylab = "Mean connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.8)
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6
  warning("No automatic soft power selected. Using 6.")
}

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

write.csv(data.frame(gene_id = colnames(datExpr), module = moduleColors),
          file.path(outdir, "gene_modules.csv"), row.names = FALSE)
write.csv(MEs, file.path(outdir, "module_eigengenes.csv"))

# convert categorical traits to numeric model matrix
trait_model <- model.matrix(~ . - 1, data = traits)
moduleTraitCor <- cor(MEs, trait_model, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

write.csv(moduleTraitCor, file.path(outdir, "module_trait_correlation.csv"))
write.csv(moduleTraitPvalue, file.path(outdir, "module_trait_pvalue.csv"))

pdf(file.path(outdir, "module_trait_heatmap.pdf"), width = 10, height = 8)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(trait_model),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

cat("WGCNA completed: ", outdir, "\n")
