# DeSEq2-Differential-MicroRNA-Expression-Analysis-of-Hepatoblastoma-GSE153089-
DeSEq2 Differential MicroRNA Expression Analysis of Hepatoblastoma (GSE153089)
# 🧬 DESeq2 DME Pipeline: MicroRNA Expression in Hepatoblastoma (GSE153089)

This R script automates a complete bioinformatics workflow for **Differential MicroRNA Expression (DME)** analysis of the **GSE153089** dataset. The study compares **Hepatoblastoma Tumor** tissue against **Nontumorous Surrounding Liver** tissue (Normal).

The pipeline uses the robust **`DESeq2`** package, adapting it for MicroRNA microarray data, and provides essential steps for automated data fetching, sample grouping, and comprehensive visualization.

## 🚀 Key Features

* **DME Analysis:** Specifically targets MicroRNA expression analysis, which often requires careful data handling.
* **DESeq2 Adaptation:** Utilizes the `DESeq2` package on rounded, normalized array values, which is a common robust workflow for DGE/DME analysis of continuous data when count-based methods are desired.
* **Automated Data Retrieval:** Downloads expression data and metadata directly from **GEO (GSE153089)** using `GEOquery`.
* **Robust Sample Grouping:** Correctly identifies and groups samples into **"Tumor"** and **"Normal"** categories based on the trailing 'N' in the sample titles.
* **Integrated Visualization:** Generates three essential plots for interpreting the results: **PCA**, **Heatmap**, and **Volcano Plot**.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE153089 | Study of a microRNA cluster in hepatoblastoma. |
| **DME Tool** | `DESeq2` (Adapted for Array Data) | Statistical method for robust differential expression testing. |
| **Comparison** | Tumor vs. Normal | Identifies microRNAs regulated in hepatoblastoma tumor tissue. |
| **Significance** | $\text{padj} < 0.05$, $\text{FCcutoff} = 1.5$ | Used for filtering significant findings. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically installs and loads the necessary Bioconductor and CRAN packages:
* `GEOquery` (For data download)
* `DESeq2` (For DME analysis)
* `pheatmap` (For Heatmap visualization)
* `EnhancedVolcano` (For Volcano Plot visualization)
* `ggplot2` (For PCA visualization)
* `matrixStats` (For efficient row-wise variance calculation)

### ⚙️ Execution

1.  **Download** the `DeSeq2 Differential MicroRNA Expression Analysis of Hepatoblastoma (GSE153089).R` file.
2.  **Optional:** The output path is set to `D:/DOWNLOADS/` by default (Step 2). You can change this path if needed.
3.  **Execute** the script in your R environment:
    ```R
    source("DeSeq2 Differential MicroRNA Expression Analysis of Hepatoblastoma (GSE153089).R")
    ```

---

## 📁 Output Files (3 Plots + 1 CSV)

All output files are saved to the specified `save_path` (default: `D:/DOWNLOADS/`).

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `Significant_DEGs_GSE153089.csv` | CSV | Table containing all microRNAs with an adjusted p-value (padj) $< 0.05$. |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `pca_plot_GSE153089.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating global clustering and separation of Tumor vs. Normal samples. |
| `heatmap_top50_GSE153089.png` | QC | **Heatmap of the Top 50 Most Variable MicroRNAs** (based on VST-transformed data) to visualize sample grouping quality. |
| `volcano_plot_GSE153089.png` | Results | **Volcano Plot** showing the $\log_2 \text{Fold Change}$ vs. $P_{\text{value}}$, highlighting significant and highly changed microRNAs. |
