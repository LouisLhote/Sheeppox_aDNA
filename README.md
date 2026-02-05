## Sheeppox_aDNA

Repository containing the analysis pipelines and scripts used in the study on ancient DNA from sheeppox virus (Capripoxvirus).  
Add your article link here: `(https://doi.org/xxx)`.

---

## Overview

This repository gathers all code used for:

- **Screening** of sequencing data for Capripoxvirus
- **Mitochondrial** analyses and competitive mapping
- **Metagenomic assembly** and scaffolding
- **Gene integrity** analyses across genomes
- **Pairwise distance** calculations and recombination-aware analyses
- **Inactivated genes** profiling
- **Phylogenetic dating with BEAST**

Each directory is a self-contained module or pipeline focused on one part of the analysis.

---

## Repository Structure

- **`Screening/`**  
  Snakemake pipeline for Capripox screening.  
  Includes:
  - **Adapter removal**
  - **Host read removal**
  - **KrakenUniq-based taxonomic screening**
  - Post-pipeline scripts to compute **E-values** for Capripox viruses.

- **`Mitochondrial/`**  
  Scripts for mitochondrial read mapping and downstream processing.  
  Includes:
  - Mapping to mitochondrial references
  - Scripts to **split competitive BAMs** into per-species BAM files.

- **`Meta_assembly/`**  
  Post-processing scripts for metagenomic assemblies.  
  Ensures that read pairs are preserved before running **SPAdes** and **RagTag**.

- **`Gene_integrity/`**  
  Pipeline for gene integrity analysis across multiple genomes.  
  Works on:
  - A set of genomes in FASTA format
  - A reference GFF annotation  
  Produces summaries of intact vs. disrupted genes.

- **`Pairwise_distance/`**  
  Scripts to compute pairwise genetic distances on Capripox alignments.  
  Features:
  - Pairwise similarity across:
    - The **full genome**
    - **Sliding windows** along the genome
  - Option to **exclude recombinant regions** via an input file
  - Ability to compute distances **between species groups** using a metadata file.

- **`Inactivated_genes_analyses/`**  
  Multiple sequence alignments of inactivated genes generated with **MACSE**.  
  Scripts to:
  - Analyze **inactivation profiles**
  - Summarize patterns of gene disruption.

- **`BEAST/`**  
  All BEAST XML input files used in the paper.  
  Includes:
  - Configured **XMLs** for different clock / tree models
  - **Aligned datasets** used as BEAST input.

---





