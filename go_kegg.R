args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)

library(org.Sc.sgd.db)
library(clusterProfiler)

# Get input arguments
input_file <- args[1]
output_file <- args[2]

# Read in differential expression results from CSV file
dea <- read.csv(input_file)
dea <- dea[dea$gene_systematic_name != 'X',]

# Convert gene IDs to Entrez IDs
entrez_ids <- bitr(dea$gene_systematic_name, fromType = "ORF",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Sc.sgd.db,
                        drop = FALSE)
dea <- merge(dea, entrez_ids, by.x = 'gene_systematic_name', by.y = 'ORF', all.x = T)

# Extract DEGs
degs <- as.data.frame(rbind(cbind(subset(subset(dea, padj <= 0.05),
                                         log2FoldChange >= 0.2), DE='up-regulated'),
                            cbind(subset(subset(dea, padj <= 0.05),
                                         log2FoldChange <= -0.2), DE='down-regulated')))


# GO/KEGG analysis
goe <- data.frame()
# kegg <- data.frame()
for (c in unique(degs$contrast)) {
  for (de in unique(degs$DE[degs$contrast == c])) {
    temp.deg <- degs[degs$contrast == c & degs$DE == de,]
    temp.goe <- enrichGO(gene          = temp.deg$ENTREZID,
                         universe      = dea$ENTREZID,
                         OrgDb         = org.Sc.sgd.db,
                         keyType       = "ENTREZID",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
    if (dim(temp.goe)[1] > 0) {
      goe <- rbind(goe, data.frame(temp.goe, contrast = c, DE = de))
    }
    
    # temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
    #                         universe     = dea$ENSEMBL,
    #                         organism     = 'sce',
    #                         pAdjustMethod = "BH",
    #                         pvalueCutoff = 0.05,
    #                         qvalueCutoff  = 0.05)
    # if (dim(temp.kegg)[1] > 0) {
    #   kegg <- rbind(kegg, data.frame(temp.kegg, contrast = c, DE = de))
    # }
  }
}

write.csv(goe, output_file)

