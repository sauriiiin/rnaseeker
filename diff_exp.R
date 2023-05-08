args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)

library(tidyverse)
library(DESeq2)

# Get input arguments
input_file <- args[1]
output_file <- args[2]
info_file <- args[3]
ref_sample <- args[4]
orfid_file <- args[5]

# Read in read counts from CSV file
sample_info <- read.csv(info_file)
sample_info$name <- factor(sample_info$sample, levels = c(ref_sample,  unique(sample_info$sample[sample_info$sample != ref_sample])))


counts_df <- read.csv(input_file)
row.names(counts_df) <- counts_df$ORF
counts_df <- counts_df[,-1]
colnames(counts_df) <- paste(sample_info$sample, sample_info$replicate, sep = '_')

keep <- rowSums(counts_df) >= 10
# keep <- rowSums(cnts < 10) == 0
counts_df <- counts_df[keep,]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts_df),
                              colData = data.frame(name = sample_info$sample,
                                                  replicate = sample_info$replicate,
                                                  sample = colnames(counts_df)),
                              design = ~ name)

# Filter low count genes and perform differential expression analysis
dds <- DESeq(dds)

# Get differential expression results table
res <- data.frame()
for (c in resultsNames(dds)[2:length(resultsNames(dds))]){
  temp <- lfcShrink(dds, coef=c, type = "normal")#, type="apeglm")
  temp$orf_name <- rownames(temp)
  rownames(temp) <- NULL
  res <- rbind(res, data.frame(temp, contrast = str_remove(c, 'name_')))
}
# res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
# res_df <- as.data.frame(res)
# res_df$orf_name <- row.names(res_df)

# Write results table to CSV file
translatome_id <- read.csv(orfid_file)
translatome_id$orf_name <- sprintf('chr%d_%d', translatome_id$chr, translatome_id$first_coord)
res <- merge(res, translatome_id[,c('orf_name','gene_systematic_name')], by = 'orf_name')

write.csv(res, file = output_file, row.names = F)
