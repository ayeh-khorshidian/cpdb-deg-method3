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
```
cpdb-deg-method3/
│
├── cpdb_deg_method3.py # Main pipeline script
├── requirements.txt # Python dependencies
├── README.md # Documentation (you are here)
└── .gitignore # Ignore venv, outputs, caches, etc.
```
---

## 📥 Requirements

Install dependencies using:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
Or install globally:
pip install -r requirements.txt
The main requirements are:
pandas
numpy
matplotlib
cellphonedb
```
---

## 📂 Expected Folder Structure Under BASE_DIR
Before running the pipeline, ensure your BASE_DIR contains:
```
BASE_DIR/
│
├── DEGs/                                # Input DEG CSV files
│     ├── <celltype1>_DEGs.csv
│     ├── <celltype2>_DEGs.csv
│     └── ...
│
├── cpdb_out/
│     ├── normalised_log_counts.h5ad     # AnnData object
│     └── metadata.tsv                   # Cell metadata
│
└── cellphonedb.zip                      # CPDB database
```
---

## 🚀 How to Run the Pipeline
1. Clone this repository
git clone https://github.com/<your-username>/cpdb-deg-method3.git
cd cpdb-deg-method3
2. Install dependencies
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
3. Edit BASE_DIR
Open cpdb_deg_method3.py and modify:
BASE_DIR = Path("/path/to/project/root").resolve()
→ Set this to your actual project root directory.
4. Run the pipeline
python cpdb_deg_method3.py
The script will automatically:
Convert DEGs
Generate meta files
Run CPDB DEG Method 3
Summarize results
Create plots
Save everything under new_DEGs_analysis_output/

---

## 🧪 Output
Under new_DEGs_analysis_output/analysis_degs_method3/ you will find:
degs_degs_method_means.txt
degs_degs_method_pvalues.txt
degs_degs_method_significant_means.txt
deg_outputs_summary.csv
significant_means_summary.csv
interaction_counts_heatmap.png
top10_bar.png
Each contrast is placed in its own folder, e.g.:
DEG_Results_2025-11-18_R+1_vs_preflight/

---

## 🧬 What Are DEGs?
DEG = Differentially Expressed Gene
A gene that shows statistically significant up- or down-regulation between two conditions.
CellPhoneDB Method 3 uses DEGs to highlight ligand–receptor interactions that are driven by transcriptional changes.

---

## 📄 License
This project is licensed under the MIT License.  
See the `LICENSE` file for details.

---

## 🙋 Questions / Improvements
If you want to:
Add advanced plots (network graphs, top ligand–receptor pairs)
Add Nextflow / Snakemake workflow version
Add Jupyter notebook examples
Add tutorial datasets
Just open an issue or message me. Happy to help refine it.
