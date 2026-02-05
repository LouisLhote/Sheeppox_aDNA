# 🐑🦠 Sheeppox_aDNA

> Repository containing the analysis pipelines and scripts used in the study on ancient DNA from sheeppox virus.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-green.svg)](https://snakemake.github.io/)


**Article**: [Link to article](https://doi.org/xxx)

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Getting Started](#-getting-started)
- [Usage](#-usage)
- [Citation](#-citation)
- [License](#-license)

---

## 🧬 Overview

This repository gathers all code used for:

- **🔍 Screening** of sequencing data for Capripoxvirus
- **🧬 Mitochondrial** analyses and competitive mapping
- **🔬 Metagenomic assembly** and scaffolding
- **🧪 Gene integrity** analyses across genomes
- **📊 Pairwise distance** calculations and recombination-aware analyses
- **🔬 Inactivated genes** profiling
- **🌳 Phylogenetic dating with BEAST**

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

## 🚀 Getting Started

### Prerequisites

- Python 3.7+
- Snakemake (for screening pipeline)
- BEAST2 (for phylogenetic analyses)
- Other dependencies as specified in each subdirectory

### Clone the repository

```bash
git clone https://github.com/<your-username>/Sheeppox_aDNA.git
cd Sheeppox_aDNA
```

### Environment / Dependencies

Each subdirectory may have its own requirements (e.g. Snakemake, BEAST, MACSE, KrakenUniq, SPAdes, RagTag, R or Python packages).  
Check the `README` or configuration files within each subfolder (e.g. `Snakefile`, `environment.yml`, `requirements.txt`, or comments in scripts) for details.

---

## 💻 Usage

- **Reproducing the screening pipeline**  
  Go to `Screening/` and follow the instructions in its documentation (or the `Snakefile` comments).

- **Running mitochondrial analyses**  
  Go to `Mitochondrial/` and use the provided mapping and BAM-splitting scripts.

- **Running BEAST analyses**  
  Go to `BEAST/`, open the XML files in BEAST, and run with the provided alignments.

---

## 📚 Citation

If you use this repository or any of the pipelines in your work, please cite:

- **Main sheeppox aDNA article**:  
  Add full citation here once available.

You should also cite the tools used (BEAST, MACSE, KrakenUniq, SPAdes, RagTag, etc.) according to their documentation.

---

## 📄 License

This project is licensed under the MIT License - see the `LICENSE` file for details.

---





