# Rich_SSN

# Single-Sample Gene Coexpression Networks in Breast Cancer Subtypes

This repository contains the R scripts and analysis pipelines used in the study of gene coexpression networks across breast cancer subtypes. The analysis focuses on single-sample coexpression networks, utilizing RNA-Seq data from *The Cancer Genome Atlas (TCGA)*. The goal of this research is to uncover structural and functional alterations in the coexpression networks of different breast cancer subtypes (Luminal A, Luminal B, Her2-enriched, Basal-like) compared to normal breast tissue.

## Study Overview
Breast cancer is characterized by significant molecular heterogeneity, which poses challenges for diagnosis and treatment. By leveraging the **LIONESS** and **ARACNe** algorithms, we construct individual-specific gene coexpression networks. These networks allow for a detailed exploration of intrachromosomal (CIS) and interchromosomal (TRANS) interactions, identifying patterns of genomic regulation and instability unique to each subtype.

## Key Objectives:
- Analyze aggregated and single-sample gene coexpression networks to uncover topological differences between breast cancer subtypes and healthy tissue.
- Assess the proportion of CIS and TRANS interactions in both aggregated and single-sample networks.
- Perform survival analysis to investigate the relationship between intrachromosomal interactions and patient outcomes.
- Identify high-degree genes and their cytoband localization to uncover potential therapeutic targets specific to each subtype.

## Repository Contents:
- **R Scripts**: All scripts used for generating figures and calculating network metrics are included in the `scripts/` directory. Each script corresponds to a specific figure or section of the article.
- **Data**: The RNA-Seq data used in this study is publicly available from the TCGA database. This repository does not include the raw data, but preprocessing steps can be found in the `data_processing/` directory.
- **Figures**: Generated figures are located in the `figures/` directory, and include visualizations of coexpression networks, Kaplan-Meier survival curves, and Jaccard index comparisons.
- **Analysis Pipelines**: Detailed pipelines for network inference (ARACNe and LIONESS), network metric calculation, and survival analysis are provided in the `analysis/` directory.

