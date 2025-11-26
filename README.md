# CellPhoneDB DEG Method 3 Pipeline

This repository contains a fully reproducible Python pipeline for running **CellPhoneDB DEG Method 3**, a ligand–receptor interaction analysis that incorporates **differentially expressed genes (DEGs)** across cell types.  
The workflow includes DEG preprocessing, metadata assembly, CellPhoneDB execution, QC summaries, and basic visualizations.

---

## 📌 Key Features

- **DEG preprocessing**
  - Converts input DEG CSV files into standardized two-column TXT files (`cell_type`, `gene`) required for CPDB Method 3.

- **Metadata assembly**
  - Generates a simple meta table linking each cell type to its DEG TXT file.

- **Automated CellPhoneDB execution**
  - Runs CellPhoneDB DEG Method 3 with user-defined paths and parameters.

- **Result summarization**
  - Automatically inspects contrast folders and summarizes:
    - Means table presence
    - P-values table presence
    - Significant interactions table presence
    - Row counts per file

- **Visualization tools**
  - Heatmap of significant interactions per contrast  
  - Top-10 contrasts bar plot  
  - Summary CSV tables

- **Portable file structure**
  - All paths are constructed using a single `BASE_DIR` to ensure portability, safety, and easy GitHub sharing.

---

## 📁 Repository Structure

