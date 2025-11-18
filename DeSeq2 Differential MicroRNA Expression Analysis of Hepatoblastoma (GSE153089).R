# Title: Differential MicroRNA Expression in Hepatoblastoma (GSE153089)
# Dataset: GSE153089 (Hepatoblastoma Tumor vs. Nontumorous Surrounding Liver)
# NOTE: This script is analyzing MicroRNA array data using DESeq2 as a DGE workflow.

# 1. Installation and Library Loading
# -------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "pheatmap", "EnhancedVolcano", "ggplot2"), update = FALSE)

# Load libraries required for the analysis
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(matrixStats) 

# 2. Define and Create Save Directory
# -------------------------------------------------
save_path <- "D:/DOWNLOADS/"
if (!dir.exists(save_path)) {
  cat("Creating directory:", save_path, "\n")
  dir.create(save_path, recursive = TRUE)
}
cat("Plots will be saved to:", save_path, "\n")

# 3. STEP 1 — Load REAL GEO DATA (GSE153089)
# -------------------------------------------------
cat("Downloading GEO dataset GSE153089... This may take a moment.\n")
gse <- getGEO("GSE153089", GSEMatrix = TRUE)
gse <- gse[[1]]

# Extract expression matrix & metadata
exprSet <- exprs(gse)
pdata <- pData(gse)

cat("Dataset dimensions (Genes/miRNAs x Samples):\n")
print(dim(exprSet))

# --- Define Sample Condition (Tumor vs Normal) ---
# The sample titles use 'N' for Nontumorous surrounding livers (Normal) 
# and lack 'N' for Tumor.
pdata$condition <- ifelse(grepl("N$", pdata$title), "Normal", "Tumor")
pdata$condition <- factor(pdata$condition, levels = c("Normal", "Tumor"))

cat("\nFinal condition table:\n")
print(table(pdata$condition))

# 4. STEP 2 — Prepare and Run DESeq2
# -------------------------------------------------
# Data is MicroRNA Array. We round the normalized values for DESeq2.
exprSet <- round(exprSet)  

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = pdata,
                              design = ~ condition)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

cat("\nRunning DESeq2 pipeline...\n")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ] # Sort results by adjusted p-value

cat("\nDESeq2 results summary:\n")
summary(res)

# 5. STEP 3 — Variance Stabilizing Transformation (VST)
# -------------------------------------------------
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# 6. STEP 4 — PCA Plot (Dimensional Reduction)
# -------------------------------------------------
cat("Generating PCA Plot...\n")
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggplot2::ggtitle("PCA: Hepatoblastoma Tumor vs Normal (GSE153089)")

# --- Save PCA Plot ---
ggsave(filename = "pca_plot_GSE153089.png", plot = pca_plot, path = save_path, width = 8, height = 6, units = "in")
cat("PCA plot saved to:", file.path(save_path, "pca_plot_GSE153089.png"), "\n")

# 7. STEP 5 — Heatmap of Top 50 Variable Genes
# -------------------------------------------------
cat("Generating Heatmap...\n")
# Selects 50 most variable genes across samples
topVarGenes <- head(order(matrixStats::rowVars(vsd_mat), decreasing = TRUE), 50)

# --- Open PNG file device to save heatmap ---
png(file.path(save_path, "heatmap_top50_GSE153089.png"), width = 800, height = 1000, res = 100)

pheatmap(vsd_mat[topVarGenes, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
         main = "Top 50 Variable Genes (GSE153089)")

# --- Close the file device ---
dev.off()
cat("Heatmap saved to:", file.path(save_path, "heatmap_top50_GSE153089.png"), "\n")

# 8. STEP 6 — Volcano Plot
# -------------------------------------------------
cat("Generating Volcano Plot...\n")

volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                pCutoff = 0.05,
                                FCcutoff = 1.5,
                                title = 'Volcano Plot: Tumor vs Normal Liver',
                                subtitle = 'Differential MicroRNA Expression',
                                legendLabels = c('NS','Log2FC','p-value','p-value & Log2FC'))

# --- Save Volcano Plot ---
ggsave(filename = "volcano_plot_GSE153089.png", plot = volcano_plot, path = save_path, width = 10, height = 8, units = "in")
cat("Volcano plot saved to:", file.path(save_path, "volcano_plot_GSE153089.png"), "\n")

# 9. STEP 7 — Save Significant Results
# -------------------------------------------------
cat("Saving significant results table...\n")
sigGenes <- as.data.frame(subset(res, padj < 0.05))
write.csv(sigGenes, file.path(save_path, "Significant_DEGs_GSE153089.csv"))

cat("\nNumber of significant genes:", nrow(sigGenes), "\n")
cat("Significant DEGs saved to:", file.path(save_path, "Significant_DEGs_GSE153089.csv"), "\n")

# --- Final Confirmation Message ---
cat("\nAnalysis complete. All files saved successfully to", save_path, "\n")
# ----------------------------------
