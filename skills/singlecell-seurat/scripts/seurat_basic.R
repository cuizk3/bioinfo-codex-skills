#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript seurat_basic.R 10x_matrix_dir outdir")
}

matrix_dir <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X(data.dir = matrix_dir)
obj <- CreateSeuratObject(counts = counts, project = "fish_gonad_scRNA", min.cells = 3, min.features = 200)

# Mitochondrial gene pattern may vary by species. Adjust pattern if needed.
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-|^MT-")

pdf(file.path(outdir, "QC_violin_before_filter.pdf"), width = 9, height = 4)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20)

pdf(file.path(outdir, "UMAP_clusters.pdf"), width = 7, height = 6)
print(DimPlot(obj, reduction = "umap", label = TRUE))
dev.off()

png(file.path(outdir, "UMAP_clusters.png"), width = 1800, height = 1500, res = 300)
print(DimPlot(obj, reduction = "umap", label = TRUE))
dev.off()

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(outdir, "cluster_markers.csv"), row.names = FALSE)

saveRDS(obj, file.path(outdir, "seurat_object.rds"))

cat("Seurat analysis completed: ", outdir, "\n")
