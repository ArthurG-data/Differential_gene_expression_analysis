# Differential Gene Expression (DGE) Analysis Comparison Report

![Language](https://img.shields.io/badge/language-R%2FPython-blue)
![Framework](https://img.shields.io/badge/framework-Bioconductor%2FedgeR%2C%20DESeq2%2C%20NOISeq%2C%20Limma%2Bvoom-orange)
![Status](https://img.shields.io/badge/status-Completed-success)

## Overview

This repository provides a comprehensive analysis of Differential Gene Expression (DGE) tools for RNA-seq data. The study compares four major pipelines—**edgeR**, **DESeq2**, **limma+voom**, and **NOISeq**—based on their ability to identify differentially expressed genes (DEGs) under various conditions. The analysis uses synthetic RNA-seq datasets to evaluate the tools' sensitivity, specificity, and computational efficiency.

---

## Technology


**RNA-Sequencing (RNA-seq):**
RNA-seq is a high-throughput sequencing method used to study gene expression. It provides a quantitative measurement of RNA molecules, offering insights into cellular biology. This study focuses on **Differential Gene Expression (DGE)**, which investigates how gene expression varies across different conditions.

**Tools Evaluated:**
1. **edgeR**  
   - Statistical Model: Negative Binomial  
   - Normalization: TMM (Trimmed Mean of M-values)  
   - Strengths: High sensitivity, effective for small sample sizes.  
   - Trade-offs: Increased false positives.  
   
2. **DESeq2**  
   - Statistical Model: Negative Binomial  
   - Normalization: Size Factors  
   - Strengths: High precision and documentation.  
   - Trade-offs: Moderate sensitivity, computationally intensive.  

3. **limma+voom**  
   - Statistical Model: Linear Models with Empirical Bayes  
   - Normalization: TMM via edgeR  
   - Strengths: Excellent runtime efficiency, low false positives.  
   - Trade-offs: Lower sensitivity for small sample sizes.  

4. **NOISeq**  
   - Statistical Model: Non-parametric Empirical  
   - Normalization: Reads Per Kilobase Million (RPKM)  
   - Strengths: Suitable for high biological variability datasets.  
   - Trade-offs: Performance declines with larger sample sizes.

---

## Methods

### Experimental Design

Synthetic RNA-seq data with predefined conditions were used:
- **Conditions:** Two groups per experiment.  
- **Samples per Condition:** 3, 6, and 9.  
- **DEG Ratios:** 1:1, 3:1, and 1:0 (upregulated to downregulated genes).  

The pipelines processed count matrices to detect DEGs. Outputs were compared based on:
- **Sensitivity**: Proportion of true positives identified.  
- **Specificity**: Ability to avoid false positives.  
- **Runtime**: Computational efficiency for varying sample sizes.

### Workflow
1. **Read Alignment:** Tools like Bowtie2 mapped raw RNA-seq reads to a reference genome.  
2. **Normalization:** Addressed biases (e.g., gene length and sequencing depth).  
3. **Statistical Testing:** Applied threshold criteria (e.g., FDR < 0.05, LogFoldChange > 0.58).  
4. **Performance Metrics:** Evaluated DEG overlap across tools and computational time.

---

## Key Findings

### Sensitivity and Specificity
- **edgeR** consistently detected the highest number of DEGs but with slightly lower specificity.
- **DESeq2** achieved a balance between sensitivity and specificity, excelling in reproducibility.
- **limma+voom** produced the fewest false positives but required larger sample sizes for robust results.
- **NOISeq** struggled with higher sample sizes due to its reliance on non-parametric methods.

<p>
  <img src="sens_vs_spec_4_models.png">
</p>

### Runtime Efficiency
- **limma+voom** was the fastest, unaffected by increasing sample sizes.
- **NOISeq** exhibited the slowest runtime for small sample sizes due to additional computational steps.

<p>
  <img src="runtime-Analysis.png">
</p>

### Concordance
- **edgeR** and **DESeq2** showed the greatest agreement, likely due to similar statistical models.
- **NOISeq** outputs were often unique and less consistent as sample sizes increased.

<p>
  <img src="9_samples_concordance.png">
</p>

---

## Rationale

The study highlights the trade-offs between sensitivity and specificity in DGE tools. While **edgeR** and **DESeq2** provide comprehensive solutions, **limma+voom** offers efficiency for large datasets, and **NOISeq** adapts to high biological variability. The findings underscore the importance of tool selection based on experimental goals and dataset characteristics.

---

## Repository Structure

```plaintext
├── data/
│   ├── Synthetic RNA-seq Count Matrices
│   ├── Metadata Files
├── results/
│   ├── DGE Results by Tool and Condition
├── figures/
│   ├── Sensitivity and Specificity Plots
│   ├── Runtime Comparisons
├── scripts/
│   ├── edgeR_analysis.R
│   ├── DESeq2_analysis.R
│   ├── limma_analysis.R
│   ├── NOISeq_analysis.R
└── README.md
