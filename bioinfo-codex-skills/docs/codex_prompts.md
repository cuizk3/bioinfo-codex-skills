# Codex 常用提示词

## 1. 初始化项目

```text
请阅读 AGENTS.md，检查本仓库结构，并告诉我目前哪些文件已经可以用于生物信息分析，哪些文件还需要补充。
```

## 2. RNA-seq 分析

```text
请使用 rnaseq skill，基于 metadata/sample_info.csv、data/raw/、data/reference/genome.fa 和 data/reference/genes.gtf 创建并运行 RNA-seq 流程。不要修改 raw data。输出结果放在 results/rnaseq。
```

## 3. 差异表达解释

```text
请读取 results/rnaseq/deseq2/DEG_*.csv，重点分析 dmrt1、foxl2、figla、cyp19a1a、amh、sox9、gsdf、wnt4、rspo1 的表达变化，并写成论文结果部分的中文描述。
```

## 4. 富集分析

```text
请使用 enrichment skill，对差异表达结果做 GO/KEGG/GSEA 分析。若没有物种 OrgDb，请使用自定义 TERM2GENE 表。
```

## 5. CRISPR sgRNA 设计

```text
请使用 crispr-sgrna skill，对 data/reference/target_gene.fasta 设计 SpCas9 sgRNA，输出前 10 个候选靶点，并解释哪个最适合做敲除。
```

## 6. WGCNA

```text
请使用 wgcna skill，用 results/rnaseq/deseq2/vst_matrix.csv 和 metadata/traits.csv 构建共表达网络，筛选与 sex 或 genotype 显著相关的模块和 hub genes。
```

## 7. 单细胞分析

```text
请使用 singlecell-seurat skill，对 data/raw/filtered_feature_bc_matrix 进行 Seurat 分析，输出 QC、UMAP、cluster markers，并结合鱼类性腺细胞 marker 进行注释。
```
