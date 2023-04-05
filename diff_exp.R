args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(DESeq2)

# Get input arguments
input_file <- args[1]
output_file <- args[2]
sample_names <- strsplit(args[3], ",")[[1]]
replicate_info <- strsplit(args[4], ",")[[1]]

# Read in read counts from CSV file
counts_df <- read_csv(input_file, col_types = cols())

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df[, sample_names],
                              colData = data.frame(condition = sample_names,
                                                  replicate = replicate_info),
                              design = ~ condition)

# Filter low count genes and perform differential expression analysis
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Get differential expression results table
res <- results(dds, contrast = c("condition", sample_names[1], sample_names[2]))
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- row.names(res_df)

# Write results table to CSV file
write_csv(res_df, output_file)
