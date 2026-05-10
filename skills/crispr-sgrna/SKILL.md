---
name: crispr-sgrna
description: Design and rank CRISPR sgRNAs for SpCas9 NGG targets from gene or CDS sequences, with GC content, position, strand, simple quality flags and fish gene-editing interpretation.
---

# CRISPR sgRNA Skill

Use this skill when the user asks for:
- sgRNA design
- CRISPR/Cas9 target screening
- target sequence selection
- mutation validation planning
- editing a fish gene such as dmrt1, figla, foxl2 or cyp19a1a

## Inputs

- FASTA sequence of target gene, exon, CDS or genomic region.
- Optional exon/CDS annotation.

## Standard output

- spacer sequence
- PAM
- strand
- start and end position
- GC content
- warnings
- rank score
- recommendation

## Rules

- For SpCas9, find N20-NGG on plus strand and CCN-N20 on minus strand.
- Prefer 40%-70% GC.
- Avoid long homopolymers.
- Avoid ambiguous bases.
- For knockout, prefer early CDS exons.
- Real off-target checking must use whole genome search. This script only performs sequence-level primary filtering.
