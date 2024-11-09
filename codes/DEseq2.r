library(DESeq2)
library(tidyr)
library(NMF)
library(brew)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# Prompt user for file paths
counts_path <- readline(prompt = "Enter the path to the counts file (e.g., 'path/to/all_mirna_counts.csv'): ")
pheno_path <- readline(prompt = "Enter the path to the phenotype data file (e.g., 'path/to/phenodata.csv'): ")
output_dir <- readline(prompt = "Enter the directory path for saving output files (e.g., 'path/to/output_dir'): ")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load the gene expression data
data <- read.csv(counts_path, row.names = 1)
dim(data)

# Load the phenotype data
pheno <- read.csv(pheno_path, row.names = 1)
table(pheno$type)

# Filter and match data and phenotype samples
meta_data_f <- pheno[which(rownames(pheno) %in% colnames(data)),]
meth_f <- data[, which(colnames(data) %in% rownames(meta_data_f))]
meta_data_f <- meta_data_f[match(colnames(meth_f), rownames(meta_data_f)),]

# Convert data to numeric matrix for further analysis
numeric_matrix <- as.matrix(data)
hist(log(numeric_matrix), col = "orange", main = "RNA (Log Scale)")

# Box plot and QQ plot for data distribution
boxplot(numeric_matrix, main = "Box Plot with Outliers", outline = TRUE)
qqnorm(data[30,])
qqline(data[30,])

# Prepare data for DESeq2 analysis
genes <- row.names(data)
data <- apply(round(data), 2, as.integer)
row.names(data) <- genes

# Reorder columns to match phenotype data
genomic_idx <- match(rownames(pheno), colnames(data))
data <- data[, genomic_idx]

# Differential Expression Analysis using DESeq2
dds <- DESeqDataSetFromMatrix(countData = data, colData = pheno, design = ~ type)
dds.run <- DESeq(dds)

# Define conditions
cond1 <- "F"
cond2 <- "UF"
cond3 <- "C"

# Perform differential expression analysis with contrasts
res1 <- results(dds.run, contrast = c("type", cond1, cond2))
res2 <- results(dds.run, contrast = c("type", cond2, cond3))
res3 <- results(dds.run, contrast = c("type", cond1, cond3))

# Filter significant results
res1 <- as.data.frame(res1[complete.cases(res1), ])
res2 <- as.data.frame(res2[complete.cases(res2), ])
res3 <- as.data.frame(res3[complete.cases(res3), ])

# Select DEGs based on adjusted p-value and log fold change
deseq.deg1 <- res1[res1$padj < 0.05,]
dgs1_up <- deseq.deg1[deseq.deg1$log2FoldChange > 1.5,]
dgs1_down <- deseq.deg1[deseq.deg1$log2FoldChange < -1.5,]

#drow DEGs volcano plot ---- F vs UF
par(mfrow=c(1,1))
with(res1, plot(log2FoldChange, -log10(padj), pch=20, main="F vs UF"))
with(subset(res1, padj<.05 & (log2FoldChange)> 1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res1, padj<.05 & (log2FoldChange)< -1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-10,y=10,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

with(res2, plot(log2FoldChange, -log10(padj), pch=20, main="UF vs C"))
with(subset(res2, padj<.05 & (log2FoldChange)> 1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res2, padj<.05 & (log2FoldChange)< -1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-3.7,y=2.5,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

with(res3, plot(log2FoldChange, -log10(padj), pch=20, main="F vs C"))
with(subset(res3, padj<.05 & (log2FoldChange)> 1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res3, padj<.05 & (log2FoldChange)< -1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=5,y=0.1,c("upregulated","downgulated"), cex=.1, bty="n", col=c("blue","red"),pch=19)
# Add a vertical line at a specific x-coordinate (e.g., x=1.5)
abline(v = 1, col = "black", lty = 2)
# Add a vertical line at a specific x-coordinate (e.g., x=1.5)
abline(v = -1, col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)

# Add a vertical line at a specific x-coordinate (e.g., x=1.5)
abline(v = 1.5, col = "black", lty = 2)
# Add a vertical line at a specific x-coordinate (e.g., x=1.5)
abline(v = -1.5, col = "black", lty = 2)


ggplot(subset_degs, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = ifelse(log2FoldChange > 1, "upregulated", "downregulated")), size = 2) +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = "F vs UF DEGs Volcano Plot", x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
    theme_minimal()

# Normalized counts and heatmap
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized = TRUE))
write.csv(as.matrix(normalized_counts), file = file.path(output_dir, "normalized_count.csv"), quote = FALSE, row.names = TRUE)

# Generate heatmap for significant DEGs
exp.degs <- as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg1), ])
pheatmap(log2(exp.degs + 1), annotation_col = meta_data_f["type"], col = rev(brewer.pal(9, "RdBu")), main = "F vs UF Heatmap", cluster_col = FALSE)


