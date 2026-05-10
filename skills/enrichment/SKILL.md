---
name: enrichment
description: Perform GO, KEGG and GSEA enrichment analysis from differential expression results, especially for fish RNA-seq and gonadal differentiation studies.
---

# Enrichment Skill

Use this skill when the user asks for:
- GO enrichment
- KEGG enrichment
- GSEA
- biological pathway interpretation
- explaining DEGs from RNA-seq

## Required inputs

- DEG table with at least:
  - gene_id
  - log2FoldChange
  - padj
- Optional gene annotation mapping:
  - gene_id to Entrez ID
  - gene_id to GO term
  - gene_id to KEGG ID

## Rules

- If the species has no built-in OrgDb, use custom TERM2GENE annotation.
- Do not over-interpret pathways from very small gene lists.
- Always report the gene universe/background used.
- For fish gonadal differentiation, highlight steroidogenesis, meiosis, germ cell development, Wnt signaling, TGF-beta signaling, retinoic acid signaling, cell cycle and apoptosis when present.
