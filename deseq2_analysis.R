# 1. Load required bioinformatics libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2", update=FALSE, ask=FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos="https://cloud.r-project.org")
}

library(DESeq2)
library(ggplot2)

print("--- Step 1: Generating Mock Multi-Sample Count Matrix ---")
set.seed(42)
genes <- paste0("gene_", 1:1000)
samples <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")

count_matrix <- matrix(
    rnbinom(6000, mu = 500, size = 10), 
    nrow = 1000, 
    dimnames = list(genes, samples)
)

# Force inject strong differentially expressed marker genes
count_matrix[1:25, 4:6] <- count_matrix[1:25, 4:6] * 5  
count_matrix[26:50, 4:6] <- count_matrix[26:50, 4:6] / 4 
count_matrix <- round(count_matrix)

print("--- Step 2: Building Experimental Design Metadata ---")
col_data <- data.frame(
    condition = factor(c(rep("control", 3), rep("treated", 3))),
    row.names = samples
)

print("--- Step 3: Launching DESeq2 Statistical Pipeline ---")
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

print("--- Step 4: Generating a Publication-Ready Volcano Plot ---")
res$significance <- "Not Significant"
res$significance[res$log2FoldChange > 1 & res$padj < 0.05] <- "Up-regulated"
res$significance[res$log2FoldChange < -1 & res$padj < 0.05] <- "Down-regulated"

volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = "blue", "Not Significant" = "gray", "Up-regulated" = "red")) +
    theme_minimal() +
    labs(
        title = "Differential Gene Expression Workspace Profiles",
        subtitle = "Target Sample Transitions (Treated vs. Control Baseline)",
        x = "Log2 Fold Change (Effect Magnitude)",
        y = "-Log10 Adjusted P-Value (Statistical Certainty)"
    ) +
    theme(legend.position = "bottom")

ggsave("volcano_plot.png", plot = volcano_plot, width = 7, height = 6, dpi = 300)
print("SUCCESS: Statistical tables calculated and 'volcano_plot.png' saved locally!")
