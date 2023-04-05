args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(tximport)

# Get input arguments
input_dir <- args[1]
output_file <- args[2]
sample_names <- strsplit(args[3], ",")[[1]]
replicate_info <- strsplit(args[4], ",")[[1]]

# Define function to read in Salmon output and generate read counts
get_read_counts <- function(salmon_file) {
  tx2gene <- read_tsv(paste0(input_dir, "/tx2gene.tsv"), col_names = c("tx_id", "gene_id"))
  txi <- tximport(files = salmon_file, type = "salmon", tx2gene = tx2gene)
  rownames(txi$counts) <- txi$genes$gene_id
  colnames(txi$counts) <- gsub(".quant", "", basename(salmon_file))
  return(txi$counts)
}

# Read in Salmon output and generate read counts for each sample
salmon_files <- list.files(path = input_dir, pattern = "quant.sf", full.names = TRUE)
counts_list <- lapply(salmon_files, get_read_counts)
counts_df <- reduce(counts_list, full_join, by = "gene_id")

# Add sample names and replicate information to counts table
counts_df$sample_name <- sample_names[match(colnames(counts_df), sample_names)]
counts_df$replicate <- replicate_info[match(colnames(counts_df), sample_names)]

# Write counts table to CSV file
write.csv(counts_df, output_file)
