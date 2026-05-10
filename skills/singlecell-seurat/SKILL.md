---
name: singlecell-seurat
description: Analyze single-cell RNA-seq data with Seurat, including QC, normalization, PCA, UMAP, clustering, marker gene discovery, and fish gonadal cell type annotation.
---

# Single-cell Seurat Skill

Use this skill when the user asks for:
- single-cell RNA-seq analysis
- Seurat workflow
- UMAP/clustering
- cell marker analysis
- gonadal cell type annotation

## Inputs

- 10x Genomics matrix directory, or
- count matrix with genes x cells

## Outputs

- QC violin plots
- UMAP
- cluster markers
- annotated Seurat object
- candidate cell type interpretation

## Fish gonad marker hints

- Germ cells: ddx4/vasa, dazl, nanos3, piwil1
- Sertoli-like cells: amh, sox9, dmrt1
- Leydig/steroidogenic cells: star, cyp11a1, hsd3b
- Granulosa-like cells: foxl2, cyp19a1a, rspo1, wnt4
- Oocyte-related cells: figla, zona pellucida genes
