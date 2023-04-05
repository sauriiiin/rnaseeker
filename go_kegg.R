args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(reactome.db)
library(org.Sc.sgd.db)
library(clusterProfiler)

# Get input arguments
input_file <- args[1]
output_file <- args[2]

# Read in differential expression results from CSV file
res_df <- read_csv(input_file, col_types = cols())

# Convert gene IDs to Entrez IDs
entrez_ids <- mapIds(org.Sc.sgd.db, keys = res_df$ensembl_gene_id,
                     keytype = "ENSEMBL", column = "ENTREZID")
res_df$entrez_id <- entrez_ids[res_df$ensembl_gene_id