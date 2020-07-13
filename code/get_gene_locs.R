library(biomaRt)

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(mart=ensembl, dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "transcription_start_site"),
               mart=ensembl)
chroms <- c(as.character(seq(1,22)))

# subset to genes on somatic chromosomes plus X (no males)
genes <- subset(genes, chromosome_name %in% chroms)
# only those that have an HGNC symbol
genes <- subset(genes, hgnc_symbol != "")
# if multiple TSS listed, go with the first
genes <- genes[!duplicated(genes$hgnc_symbol),]

saveRDS(genes, "../data/biomart_genepos.rds")
