rm(list = ls())

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(Matrix)
library(slingshot)

# this script preprocesses the single cell expression data into pseudobulk 
# either for all cells or for cells of a single lineage assigned by slingshot
# formats it into a data.table, where column values are:
# [genes, expression value, cell line, day]
# produces expression PCs and cell line PCs based on the preprocessed data

## set hyperparameters
# we can choose what kind of data we want to extract
# if `all`, then pseudobulk is for all cells
# if `i` where `i` is a string, then pseudobulk is for cells belonging to the ith lineage
which_data <- "2"
# subset samples (a cell line/day combo) to those with at least `min_reads` number of reads 
# for now we are using 50K for all cells, 10K for lineage 2
min_reads <- 10000

## read in data
# read in seurat data for a subset of 20K cells
sc <- readRDS("data/mitofilt_singlets.pca.cc.pseudo.rds")
# get snp and genotype matrices
snps <- read.table("data/genotypes_05cut.txt", header=T, row.names=1)
snp_locs <- read.table("data/snp_locations_05cut.txt", header=T, row.names=1)
# slingshot results on this subset of cells
ss <- readRDS("data/slingshot.rds")
# set to the number of lineages return by slingshot
num_lineages <- 3

## subset the data to only cells belonging to the ith lineage
if (which_data != "all") {
  ## get cells of the ith lineage
  branch_id <- slingBranchID(ss)
  
  # check if `i` is an element of a branch id string `id_str`
  # e.g. has_branch_i("1,2", "1") == TRUE, has_branch_i("1,2", 3) == FALSE
  has_branch_i <- function(id_str, i) { 
    ids <- unlist(strsplit(id_str, ","))
    return(i %in% ids)
  }
  
  ids_with_branch <- levels(branch_id)[sapply(levels(branch_id), has_branch_i, i=which_data)]
  # the subset of cells in the ith lineage
  cells_subset <- names(branch_id[branch_id %in% ids_with_branch])
  sc <- sc[, cells_subset]
}

## get raw expression data and metadata
counts <- t(sc@assays[["RNA"]]@counts)
days <- sc$diffday
inds <- factor(sc$demux.sng1)
# use cell line id and day to create unique id for each sample
sample_ids <- c()
for (i in levels(inds)) {
  for (d in levels(days)) {
    sample_ids <- append(sample_ids, paste0(i, '_', strsplit(d, ' ')[[1]][2]))
  }
}

## subset counts to what appears in biomart
# this section can be removed after full seurat obj use
hgnc <- readRDS("data/biomart_genepos.rds")
hgnc <- subset(hgnc, hgnc_symbol %in% colnames(counts))
counts <- counts[,colnames(counts) %in% hgnc$hgnc_symbol]

## get pseudobulk for each ind/day sample
tables <- c()
for (id in levels(inds)) {
  for (day in levels(days)) {
    cells <- rownames(inds)[inds == id & days == day]
    # otherwise can't do colSums
    if (length(cells) < 2) {next}
    expr <- colSums(counts[cells,])
    pseudobulk <- data.table(id = id, day = day, genes = names(expr), value = expr)
    tables <- append(tables, list(pseudobulk))
  }
}
raw_pseudobulk <- rbindlist(tables)
raw_pseudobulk$sample_id <- apply(raw_pseudobulk, 1, function(row){paste0(row[1], '_', strsplit(row[2], ' ')[[1]][2])})
saveRDS(raw_pseudobulk, paste0("data/dynamic_pb/lineage", which_data, "_raw_pseudobulk.rds"))
rm(tables, pseudobulk, cells, expr)

## subset to samples with at least certain number of reads
reads <- raw_pseudobulk %>% group_by(id, day) %>% select(sample_id, value) %>% summarise(sample_id = first(sample_id), sum=sum(value))
# plot the number of reads for each sample
p <- ggplot(data = reads, aes(x=sample_id, y=sum)) + geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot = p, width = 32, height = 8, filename = paste0("figs/dynamic_pb/lineage", which_data, "_pseudobulk_reads.png"))
# sample ids that passed the filter
sample_id_f <- reads$sample_id[reads$sum>=min_reads]
# sample ids that didn't pass the filter
sample_id_nf <- setdiff(sample_ids, sample_id_f)
raw_pseudobulk <- raw_pseudobulk[sample_id %in% sample_id_f]

## subset to snps with MAF>=0.05
afs <- rowSums(snps, na.rm=TRUE)/(2*ncol(snps))
mafs <- sapply(afs, function(x){min(x,1-x)})
mafs <- mafs[mafs>=0.05]
snps_f <- snps[names(mafs),]
snp_locs_f <- snp_locs[names(mafs),]

## subset to genes with nonzero variance
genes_sd <- raw_pseudobulk %>% group_by(genes) %>% select(sample_id, value) %>% summarise(sd = sd(value))
raw_pseudobulk <- raw_pseudobulk[genes %in% genes_sd$genes]

## subset to genes with at least 5 reads in at least 5 samples
# the number of samples in which each gene has at least 5 reads
nreads5_genes <- raw_pseudobulk[, lapply(.SD, function(x){sum(x>=5)}), by = genes, .SDcols = 'value']
# the genes with at least 5 reads in at least 5 samples
nsamples5_genes <- nreads5_genes[value >= 5]
raw_pseudobulk <- raw_pseudobulk[genes %in% nsamples5_genes$genes]

## log normalize within each sample, across all genes
pseudobulk_lognorm <- raw_pseudobulk[, lapply(.SD, function(x){log(10^4*x/sum(x)+1)}), by = sample_id, .SDcols = 'value']
pseudobulk_lognorm <- pseudobulk_lognorm[, c("id", "day", "gene") := list(raw_pseudobulk$id, raw_pseudobulk$day, raw_pseudobulk$genes)]
saveRDS(pseudobulk_lognorm, paste0("data/dynamic_pb/lineage", which_data, "_lognorm_pseudobulk.rds"))

## center and scale within each gene, across all samples
pseudobulk_standardized <- pseudobulk_lognorm[, lapply(.SD, function(x){x-mean(x, na.rm = TRUE)/sd(x, na.rm = TRUE)}), by = gene, .SDcols = 'value']
index <- pseudobulk_lognorm[, .SD, by = gene, .SDcols = c('sample_id', 'id', 'day')]
pseudobulk_standardized <- pseudobulk_standardized[, c('sample_id', 'id', 'day') := list(index$sample_id, index$id, index$day)]
rm(index)
pseudobulk_standardized <- pseudobulk_standardized[order(id, day, gene)]
saveRDS(pseudobulk_standardized, paste0("data/dynamic_pb/lineage", which_data, "_standardized_pseudobulk.rds"))

# `pseudobulk` is the preprocessed data after all the steps except for standardization
# right now the last step is log normalization
# since different standardization is needed for expression PCs and cell line PCs
pseudobulk <- pseudobulk_lognorm

## get expression PCs
# wide format for pca, note the first column is sample_id, since data.table doesn't support rownames
pseudobulk_pca_wf <- dcast(pseudobulk, sample_id + id + day ~ gene, value.var = "value")
pca_index <- pseudobulk_pca_wf[, c('sample_id', 'id', 'day')]
pseudobulk_pca_wf <- scale(pseudobulk_pca_wf[, !c("sample_id", "id", "day")])
expr_usv <- svd(pseudobulk_pca_wf)
pc_percent <- ((expr_usv$d^2)/sum((expr_usv$d^2)))*100
# sanity check variance explained and pc1 vs pc2
plot(pc_percent, pch = 20, xlim = c(1,50), xlab = "PC Number", 
     ylab = "% Variance explained", main = "Scree plot")
pca_index[, c("PC1", "PC2") := list(expr_usv$u[, 1], expr_usv$u[, 2])]
p <- ggplot(pca_index, aes(x = PC1, y = PC2, color = day)) + geom_point()
ggsave(plot = p, filename = paste0("figs/dynamic_pb/pc1vs2_lineage", which_data, ".png"))
# get the first 50 PCs
expr_pc <- expr_usv$u[, 1:50]
saveRDS(expr_pc, paste0("data/dynamic_pb/lineage", which_data, "_exprPC50.rds"))

## get cell line PCs
# here we need to use the rearranged expression matrix to center and scale
# first we need to drop the cell lines for which we don't have complete samples
missing_celllines <- unique(sapply(sample_id_nf, function(x){strsplit(x, '_')}[[1]][1]))
pseudobulk_filtered <- pseudobulk[!id %in% missing_celllines,]

# drop the genes that have zero variance for all the remaining cell lines for any day, 
# otherwise we can't standardize that column
gene_day_sd <- pseudobulk_filtered %>% group_by(gene, day) %>% summarise(sd = sd(value))
gene_day_sd <- gene_day_sd %>% group_by(gene) %>% filter(all(sd > 0))
genes_clpca <- unique(gene_day_sd$gene)
pseudobulk_filtered <- pseudobulk_filtered[gene %in% genes_clpca,]

# wide format for cell line PCA
pseudobulk_clpca_wf <- dcast(pseudobulk_filtered, id ~ day + gene, value.var = "value")
clpca_index <- pseudobulk_clpca_wf[, "id"]
pseudobulk_clpca_wf <- scale(pseudobulk_clpca_wf[, !"id"])
cellline_usv <- svd(pseudobulk_clpca_wf)
clpc_percent <- ((cellline_usv$d^2)/sum((cellline_usv$d^2)))*100
plot(clpc_percent, pch = 20, xlim = c(1,13), xlab = "cell line PC Number", 
     ylab = "% Variance explained", main = "Scree plot")
cellline_pc <- cellline_usv$u[, 1:10]
saveRDS(cellline_pc, paste0("data/dynamic_pb/lineage", which_data, "_celllinePC10.rds"))


