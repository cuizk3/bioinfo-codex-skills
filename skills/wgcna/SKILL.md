---
name: wgcna
description: Run WGCNA co-expression network analysis from normalized expression matrices and trait metadata, identifying modules, module-trait correlations, hub genes and fish gonadal differentiation interpretation.
---

# WGCNA Skill

Use this skill when the user asks for:
- WGCNA
- co-expression network
- module-trait relationship
- hub genes
- gene regulatory network from RNA-seq

## Inputs

- normalized expression matrix, genes as rows and samples as columns
- trait metadata, samples as rows
- optional target trait such as sex, genotype, stage, gonad type

## Rules

- Filter low-variance genes before WGCNA.
- Check sample clustering before network construction.
- Choose soft threshold using scale-free topology.
- Interpret modules cautiously: WGCNA shows co-expression, not causality.
