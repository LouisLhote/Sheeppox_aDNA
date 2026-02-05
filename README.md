# 🐑🦠 Sheeppox_aDNA

> Repository containing the analysis pipelines and scripts used in **3,500 years of sheeppox virus evolution inferred from archaeological and codicological genomes.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-green.svg)](https://snakemake.github.io/)


**Article**: [Link to article](https://doi.org/xxx)

---

## 🧬 Overview

This repository gathers all code used for:

- **Screening** of sequencing data for Capripoxvirus
- **Mitochondrial** analyses and competitive mapping
- **SPPV_ML** Phylogenetic analysis of the aSPPV dataset with Pathphynder placement
- **Metagenomic assembly** and scaffolding
- **Gene integrity** analyses across genomes
- **Pairwise distance** calculations and recombination-aware analyses
- **Inactivated genes** profiling
- **Phylogenetic dating with BEAST**

Each directory is a self-contained module or pipeline focused on one part of the analysis.

---

## 📁 Repository Structure

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

- **`SPPV_ML/`**  
  Phylogenetic analysis of the aSPPV dataset with **PathPhynder** placement.  
  Includes:
  - Multiple sequence alignments for SPPV genomes
  - ML trees
  - PathPhynder outputs placement

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

