"""
CellPhoneDB DEG(Differentially Expressed Gene) Method 3 pipeline (script version)

This script:
1) Converts DEG CSV files to two-column TXT files for CellPhoneDB.
2) Builds a meta file for DEG TXT files.
3) Runs CellPhoneDB DEG Method 3.
4) Summarizes outputs and generates simple plots.

All paths are defined relative to BASE_DIR so the code is portable.
"""

# =============================================================================
# 0. Imports & Global Configuration
# =============================================================================

from pathlib import Path
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cellphonedb.src.core.methods import cpdb_degs_analysis_method


# -----------------------------------------------------------------------------
# Base project directory (EDIT THIS FOR YOUR ENVIRONMENT)
# -----------------------------------------------------------------------------
BASE_DIR = Path("/path/to/project/root").resolve()

# =============================================================================
# 1. Paths & Folder Structure
# =============================================================================

# Input paths
DEG_CSV_DIR = BASE_DIR / "DEGs"
DEG_TXT_DIR = DEG_CSV_DIR / "I_made_these_txt_files"
DEG_META_FILE = DEG_TXT_DIR / "meta_method3.txt"

COUNTS_FILE = BASE_DIR / "cpdb_out" / "normalised_log_counts.h5ad"
META_FILE = BASE_DIR / "cpdb_out" / "metadata.tsv"
CPDB_ZIP = BASE_DIR / "cellphonedb.zip"

# Output paths
OUTPUT_ROOT = BASE_DIR / "new_DEGs_analysis_output"
DEG_OUTPUT_DIR = OUTPUT_ROOT / "analysis_degs_method3"


# =============================================================================
# 2. Step 1 – Convert DEG CSV files to 2-column TXT
# =============================================================================

def convert_deg_csv_to_txt():
    """
    Convert all *.csv DEG files in DEG_CSV_DIR to two-column TXT files
    (cell_type, gene) in DEG_TXT_DIR for CellPhoneDB Method 3.
    """
    DEG_TXT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[STEP 1] Converting CSV DEGs in {DEG_CSV_DIR} to TXT in {DEG_TXT_DIR}")

    for csv_path in sorted(DEG_CSV_DIR.glob("*.csv")):
        print(f"  - Processing: {csv_path.name}")

        df = pd.read_csv(csv_path)

        # Assume the first column is the gene column
        gene_col = df.columns[0]

        # Extract cell type from file name, e.g. "CD16Mono_DEGs.csv" -> "CD16Mono"
        fname = csv_path.name
        cell_type = fname.replace("_DEGs.csv", "").replace(".csv", "")

        df["cell_type"] = cell_type

        # Keep only the two columns required for CellPhoneDB
        txt_df = df[["cell_type", gene_col]].dropna().drop_duplicates()

        out_name = fname.replace(".csv", "_DEGs.txt")
        out_path = DEG_TXT_DIR / out_name

        # Save as tab-separated, without header
        txt_df.to_csv(out_path, sep="\t", index=False, header=False)

        print(f"    -> Saved: {out_path.name}")

    print("[STEP 1] Done.\n")


# =============================================================================
# 3. Step 2 – Build meta file for DEG TXT files (meta_method3.txt)
# =============================================================================

def build_deg_meta_file():
    """
    Build a simple meta file with:
        celltype, deg_txt_path
    This is useful for tracking which TXT file belongs to which cell type.
    """
    DEG_TXT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[STEP 2] Building DEG meta file at {DEG_META_FILE}")

    celltypes = []
    filepaths = []

    for txt_file in sorted(DEG_TXT_DIR.glob("*_DEGs.txt")):
        name = txt_file.name
        celltype = name.split("_DEGs")[0]
        celltypes.append(celltype)
        filepaths.append(str(txt_file.resolve()))

    df_meta = pd.DataFrame({"celltype": celltypes, "deg_txt_path": filepaths})
    df_meta.to_csv(DEG_META_FILE, sep="\t", index=False)

    print(df_meta.head())
    print(f"[STEP 2] Saved meta file with {len(df_meta)} entries.\n")


# =============================================================================
# 4. Step 3 – Configure and Run CellPhoneDB DEG Method 3
# =============================================================================

def run_cpdb_degs_method3(
    counts_file: Path = COUNTS_FILE,
    meta_file: Path = META_FILE,
    cpdb_zip: Path = CPDB_ZIP,
    deg_txt_dir: Path = DEG_TXT_DIR,
    output_dir: Path = DEG_OUTPUT_DIR,
    counts_data: str = "log-normalized",
    threads: int = 8,
):
    """
    Run cpdb_degs_analysis_method for DEG Method 3.

    Parameters
    ----------
    counts_file : Path
        Normalized AnnData (.h5ad) file (genes x cells).
    meta_file : Path
        Metadata TSV file with per-cell annotations.
    cpdb_zip : Path
        CellPhoneDB database zip file.
    deg_txt_dir : Path
        Folder containing two-column TXT DEG files.
    output_dir : Path
        Folder where DEG outputs will be written.
    counts_data : str
        Type of counts (e.g., 'log-normalized').
    threads : int
        Number of threads for CellPhoneDB.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    print("[STEP 3] Running CellPhoneDB DEG Method 3 with the following config:")
    print(f"  COUNTS_FILE : {counts_file}")
    print(f"  META_FILE   : {meta_file}")
    print(f"  CPDB_ZIP    : {cpdb_zip}")
    print(f"  DEG_TXT_DIR : {deg_txt_dir}")
    print(f"  OUTPUT_DIR  : {output_dir}")
    print()

    try:
        cpdb_degs_analysis_method.call(
            meta_file_path=str(meta_file),
            counts_file_path=str(counts_file),
            database_file_path=str(cpdb_zip),
            degs_folder_path=str(deg_txt_dir),
            counts_data=counts_data,
            log_level="INFO",
            threads=threads,
            output_path=str(output_dir),
        )
        print("[STEP 3] CellPhoneDB DEG analysis completed successfully.\n")

    except Exception as e:
        print("[STEP 3] ERROR while running CellPhoneDB DEG analysis:")
        print(e)
        print()
        raise


# =============================================================================
# 5. Step 4 – Summarize DEG outputs per contrast
# =============================================================================

def summarize_cpdb_degs_outputs(deg_output_root: Path = DEG_OUTPUT_DIR) -> pd.DataFrame:
    """
    For each contrast folder (e.g. DEG_Results_YYYY-MM-DD_R+1_vs_preflight),
    check the presence of:
        - degs_degs_method_means.txt
        - degs_degs_method_pvalues.txt
        - degs_degs_method_significant_means.txt

    Returns
    -------
    summary_df : pd.DataFrame
        One row per contrast with basic file existence and row counts.
    """
    deg_output_root = Path(deg_output_root)
    print(f"[STEP 4] Summarizing DEG outputs under: {deg_output_root}")

    summary_records = []

    for contrast_dir in sorted(deg_output_root.iterdir()):
        if not contrast_dir.is_dir():
            continue

        contrast_name = contrast_dir.name
        print(f"  - Contrast: {contrast_name}")

        means_path = contrast_dir / "degs_degs_method_means.txt"
        pvalues_path = contrast_dir / "degs_degs_method_pvalues.txt"
        sig_means_path = contrast_dir / "degs_degs_method_significant_means.txt"

        files_exist = {
            "means": means_path.exists(),
            "pvalues": pvalues_path.exists(),
            "sig_means": sig_means_path.exists(),
        }
        print(f"      Files present: {files_exist}")

        record = {
            "contrast": contrast_name,
            "means_exists": files_exist["means"],
            "pvalues_exists": files_exist["pvalues"],
            "sig_means_exists": files_exist["sig_means"],
            "means_rows": None,
            "pvalues_rows": None,
            "sig_means_rows": None,
        }

        if files_exist["means"]:
            means_df = pd.read_csv(means_path, sep="\t")
            record["means_rows"] = len(means_df)

        if files_exist["pvalues"]:
            pvalues_df = pd.read_csv(pvalues_path, sep="\t")
            record["pvalues_rows"] = len(pvalues_df)

        if files_exist["sig_means"]:
            sig_means_df = pd.read_csv(sig_means_path, sep="\t")
            record["sig_means_rows"] = len(sig_means_df)

        summary_records.append(record)

    summary_df = pd.DataFrame(summary_records)
    print("\n[STEP 4] Summary of contrast-level DEG outputs:")
    print(summary_df)
    print()

    summary_path = deg_output_root / "deg_outputs_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[STEP 4] Saved summary table to: {summary_path}\n")

    return summary_df


# =============================================================================
# 6. Step 5 – Summarize "significant_means" & Simple Plots
# =============================================================================

def summarize_cpdb_significant_means(deg_output_root: Path = DEG_OUTPUT_DIR) -> pd.DataFrame:
    """
    For each contrast, load 'degs_degs_method_significant_means.txt' and compute
    the number of significant interactions.

    Returns
    -------
    summary_df : pd.DataFrame
        One row per contrast with n_rows and total_significant_flags.
    """
    deg_output_root = Path(deg_output_root)
    print(f"[STEP 5] Summarizing significant means under: {deg_output_root}")

    records = []

    for contrast_dir in sorted(deg_output_root.iterdir()):
        if not contrast_dir.is_dir():
            continue

        contrast_name = contrast_dir.name
        sig_means_path = contrast_dir / "degs_degs_method_significant_means.txt"

        if not sig_means_path.exists():
            print(f"  [WARNING] Missing sig_means file in {contrast_name}")
            continue

        df = pd.read_csv(sig_means_path, sep="\t")

        sig_cols = [c for c in df.columns if c.startswith("significant")]
        sig_counts = df[sig_cols].sum()
        total_sig = sig_counts.sum()

        records.append({
            "contrast": contrast_name,
            "n_rows": len(df),
            "total_significant_flags": total_sig,
        })

    summary_df = pd.DataFrame(records)
    print("\n[STEP 5] Summary of significant means per contrast:")
    print(summary_df)
    print()

    out_path = deg_output_root / "significant_means_summary.csv"
    summary_df.to_csv(out_path, index=False)
    print(f"[STEP 5] Saved significant means summary to: {out_path}\n")

    return summary_df


def plot_interaction_count_heatmap(summary_df: pd.DataFrame, out_dir: Path = DEG_OUTPUT_DIR):
    """
    Create a simple heatmap of the number of significant interactions per contrast.
    """
    if summary_df.empty:
        print("[PLOT] Empty summary_df, skipping heatmap.")
        return

    out_dir = Path(out_dir)

    data = summary_df.set_index("contrast")[["n_rows"]]

    plt.figure(figsize=(6, 4))
    plt.imshow(data.values, aspect="auto")
    plt.colorbar(label="Number of significant interactions")
    plt.yticks(range(len(data.index)), data.index)
    plt.xticks([0], ["n_rows"])
    plt.title("Number of significant interactions per contrast")

    plt.tight_layout()
    out_path = out_dir / "interaction_counts_heatmap.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"[PLOT] Saved heatmap to: {out_path}")


def plot_top10_bar(summary_df: pd.DataFrame, out_dir: Path = DEG_OUTPUT_DIR):
    """
    Create a bar plot of the top 10 contrasts by number of significant interactions.
    """
    if summary_df.empty:
        print("[PLOT] Empty summary_df, skipping bar plot.")
        return

    out_dir = Path(out_dir)

    top10 = summary_df.sort_values("n_rows", ascending=False).head(10)

    plt.figure(figsize=(8, 4))
    plt.bar(top10["contrast"], top10["n_rows"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of significant interactions")
    plt.title("Top 10 contrasts by significant interactions")

    plt.tight_layout()
    out_path = out_dir / "top10_bar.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"[PLOT] Saved bar plot to: {out_path}")


# =============================================================================
# 7. Main entry point
# =============================================================================

def main():
    """
    Run the full pipeline:
    1) Build TXT DEG files.
    2) Build meta_method3 for DEGs.
    3) Run CellPhoneDB DEG Method 3.
    4) Summarize outputs.
    5) Summarize significant means and generate plots.
    """
    print("=== CellPhoneDB DEG Method 3 Pipeline ===")
    print(f"BASE_DIR: {BASE_DIR}\n")

    convert_deg_csv_to_txt()
    build_deg_meta_file()
    run_cpdb_degs_method3()
    summarize_cpdb_degs_outputs()
    sig_summary = summarize_cpdb_significant_means()
    plot_interaction_count_heatmap(sig_summary)
    plot_top10_bar(sig_summary)

    print("\n=== Pipeline finished. ===")


if __name__ == "__main__":
    main()
  
