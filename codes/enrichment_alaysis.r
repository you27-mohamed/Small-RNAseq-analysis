# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GOplot", quietly = TRUE)) BiocManager::install('GOplot')
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichR", quietly = TRUE)) BiocManager::install('enrichR')

# Load libraries
library("enrichR")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("GOplot")

# Set working directory and file paths based on user input
working_dir <- readline(prompt = "Enter the working directory path: ")
setwd(working_dir)

data_file <- readline(prompt = "Enter the path to the gene data file (e.g., 'path/to/downregulated_genes.txt'): ")
output_dir <- readline(prompt = "Enter the output directory for saving results (e.g., 'path/to/output_dir'): ")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Define the databases for enrichment analysis
dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
dbs_pw <- c("KEGG_2021_Human", "WikiPathway_2023_Human")
dbs_dd <- c("PheWeb_2019", "ClinVar_2019")

# Perform enrichment analysis
data <- read.csv(data_file)
upEnriched_go <- enrichr(genes = data$Gene.ID, databases = dbs_go)
upEnriched_pw <- enrichr(genes = data$Gene.ID, databases = dbs_pw)

# Save enrichment results to the specified output directory
write.csv(upEnriched_go$GO_Molecular_Function_2023, file.path(output_dir, "enrichment_MF_PIWI_downregulated.csv"))
write.csv(upEnriched_go$GO_Cellular_Component_2023, file.path(output_dir, "enrichment_CC_PIWI_downregulated.csv"))
write.csv(upEnriched_go$GO_Biological_Process_2023, file.path(output_dir, "enrichment_BP_PIWI_downregulated.csv"))
write.csv(upEnriched_pw$KEGG_2021_Human, file.path(output_dir, "enrichment_KEGG_PIWI_downregulated.csv"))
write.csv(upEnriched_pw$WikiPathway_2023_Human, file.path(output_dir, "enrichment_WP_PIWI_downregulated.csv"))

# Define a function to plot significant enriched terms with filtering
plotSignificantEnrich <- function(enrich_data, title, showTerms, numChar, y, orderBy) {
  significant_terms <- enrich_data[enrich_data$P.value < 0.05, ]
  
  plotEnrich(significant_terms, showTerms = showTerms, numChar = numChar, y = y, orderBy = orderBy) +
    ggtitle(title)
}

# Generate enrichment plots and save them in the specified directory
png(file.path(output_dir, "Molecular_Function_Enrichment.png"), width = 10, height = 8, units = "in", res = 300)
plotSignificantEnrich(upEnriched_go[[1]], title = "Molecular Function Enrichment", showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value")
dev.off()

svg(file.path(output_dir, "Cellular_Component_Enrichment.svg"), width = 10, height = 8)
plotSignificantEnrich(upEnriched_go[[2]], title = "Cellular Component Enrichment", showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value")
dev.off()

svg(file.path(output_dir, "Biological_Process_Enrichment.svg"), width = 10, height = 8)
plotSignificantEnrich(upEnriched_go[[3]], title = "Biological Process Enrichment", showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value")
dev.off()

svg(file.path(output_dir, "KEGG_Pathway_Enrichment.svg"), width = 10, height = 8)
plotSignificantEnrich(upEnriched_pw[[1]], title = "KEGG Pathway Enrichment", showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value")
dev.off()

svg(file.path(output_dir, "WikiPathways_Enrichment.svg"), width = 10, height = 8)
plotSignificantEnrich(upEnriched_pw[[2]], title = "WikiPathways Enrichment", showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value")
dev.off()
