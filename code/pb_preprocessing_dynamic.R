library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(Matrix)
library(slingshot)
library(preprocessCore)

# This script preprocesses the single cell expression data into pseudobulk 
# produces expression PCs and cell line PCs based on the preprocessed data
# pseudobulk aggregation can be for: 1. all cells; 2. for cells of a single lineage assigned by slingshot
# time can be: 1. collection day; 2. pseudotime deciles (only for lineage specific cells)
# formats it into a data.table, where column values are:
# [gene, expression value, cell line, time]

###############################
## set hyperparameters
###############################

# if "all", then pseudobulk is for all cells
# if `i`` where i` is a string, then pseudobulk is for cells belonging to the ith lineage
which_data <- "all"
# if "day", then use collection day as time
# if "pseudotime", then use pseudotime deciles as time (only for lineage specific cells)
which_time <- "day"
# subset samples (a cell line/time combo) to those with at least `min_reads` number of reads 
min_reads <- 50000
# directories for output
data_dir <- "data/dynamic_pb/"
fig_dir <- "figs/dynamic_pb/"
# number of expression PCs to keep
num_expr_pc <- 30
# number of cell line PCs to keep
num_cellline_pc <- 10

###############################
## read input data
###############################

# read in seurat data for a subset of 20K cells
sc <- readRDS("data/mitofilt_singlets.pca.cc.pseudo.rds")
# get snp and genotype matrices
snps <- read.table("data/genotypes_05cut.txt", header=T, row.names=1)
snp_locs <- read.table("data/snp_locations_05cut.txt", header=T, row.names=1)
# slingshot results on this subset of cells
# optional, only if we need lineage specific cells 
ss <- readRDS("data/slingshot.rds")
# set to the number of lineages return by slingshot
num_lineages <- 3

###############################
## process input data
###############################

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
  # ss <- ss[, cells_subset]
  # map the pseudotime of each cell to a decile
  pseudotime_deciles <- ntile(slingPseudotime(ss)[cells_subset, paste0("curve", which_data)], 10)
}

## get raw expression data and metadata
counts <- t(sc@assays[["RNA"]]@counts)
days <- sc$diffday
inds <- factor(sc$demux.sng1)
if (which_time == "day") {
  day_map <- c("Day 0"=0, "Day 1"=1, "Day 3"=3, "Day 5"=5, "Day 7"=7, "Day 11"=11, "Day 15"=15)
  time <- day_map[days]
  names(time) <- names(days)
} else if (which_time == "pseudotime") {
  time <- pseudotime_deciles
}
# use cell line id and time to create unique id for each sample
sample_ids <- c()
for (i in levels(inds)) {
  for (t in sort(unique(time))) {
    sample_ids <- append(sample_ids, paste0(i, '_', t))
  }
}

## subset counts to what appears in biomart
# this section can be removed after full seurat obj use
hgnc <- readRDS("data/biomart_genepos.rds")
hgnc <- subset(hgnc, hgnc_symbol %in% colnames(counts))
counts <- counts[,colnames(counts) %in% hgnc$hgnc_symbol]

## subset to snps with MAF>=0.05
afs <- rowSums(snps, na.rm=TRUE)/(2*ncol(snps))
mafs <- sapply(afs, function(x){min(x,1-x)})
mafs <- mafs[mafs>=0.05]
snps_f <- snps[names(mafs),]
snp_locs_f <- snp_locs[names(mafs),]

###############################
## define preprocessing functions 
###############################

## get raw pseudobulk from expression data
get_raw_pseudobulk <- function() {
  # get pseudobulk for each ind/time sample
  tables <- c()
  for (id in levels(inds)) {
    for (t in sort(unique(time))) {
      cells <- names(inds)[inds == id & time == t]
      # otherwise can't do colSums
      if (length(cells) < 2) {next}
      expr <- colSums(counts[cells,])
      pseudobulk <- data.table(id = id, time = t, sample_id = paste0(id, "_", t), gene = names(expr), value = expr)
      tables <- append(tables, list(pseudobulk))
    }
  }
  raw_pseudobulk <- rbindlist(tables)
  return(raw_pseudobulk)
}

## subset to samples with at least certain number of reads
subset_samples <- function(pseudobulk) {
  reads <- pseudobulk %>% group_by(id, time) %>% select(sample_id, value) %>% summarise(sample_id = first(sample_id), sum=sum(value))
  # plot the number of reads for each sample
  p <- ggplot(data = reads, aes(x=sample_id, y=sum)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plot = p, width = 32, height = 8, filename = paste0(fig_dir, "lineage", which_data, "_pseudobulk_reads.png"))
  # sample ids that passed the filter
  sample_id_f <- reads$sample_id[reads$sum>=min_reads]
  # sample ids that didn't pass the filter
  sample_id_nf <- setdiff(sample_ids, sample_id_f)
  pseudobulk <- pseudobulk[sample_id %in% sample_id_f]
  return(list(pb = pseudobulk, missing_samples = sample_id_nf))
}

## subset to genes with nonzero variance and with at least 5 reads in at least 5 samples
subset_genes <- function(pseudobulk) {
  genes_sd <- pseudobulk %>% group_by(genes) %>% select(sample_id, value) %>% summarise(sd = sd(value))
  pseudobulk <- pseudobulk[genes %in% genes_sd$genes]
  
  # the number of samples in which each gene has at least 5 reads
  nreads5_genes <- pseudobulk[, lapply(.SD, function(x){sum(x>=5)}), by = genes, .SDcols = 'value']
  # the genes with at least 5 reads in at least 5 samples
  nsamples5_genes <- nreads5_genes[value >= 5]
  pseudobulk <- pseudobulk[genes %in% nsamples5_genes$genes]
  return(pseudobulk)
}

## log normalize within each sample, across all genes
log_normalize <- function(pseudobulk) {
  pseudobulk_lognorm <- pseudobulk[, lapply(.SD, function(x){log(10^4*x/sum(x)+1)}), by = sample_id, .SDcols = 'value']
  pseudobulk_lognorm <- pseudobulk_lognorm[, c("id", "time", "gene") := list(pseudobulk$id, pseudobulk$time, pseudobulk$gene)]
  return(pseudobulk_lognorm)
}

## quantile normalize to make the gene expression distribution of each sample identical in statistical properties
quantile_normalize <- function(pseudobulk) {
  pseudobulk_wf <- dcast(pseudobulk, sample_id + id + time ~ gene, value.var = "value")
  index <- pseudobulk_wf[, c('sample_id', 'id', 'time')]
  gene_names <- colnames(pseudobulk_wf)[-c(1:3)]
  pseudobulk_qtnorm <- t(normalize.quantiles(t(pseudobulk_wf[, !c("sample_id", "id", "time")])))
  colnames(pseudobulk_qtnorm) <- gene_names
  pseudobulk_qtnorm <- data.table(pseudobulk_qtnorm)
  pseudobulk_qtnorm[, c("sample_id", "id", "time") := list(index$sample_id, index$id, index$time)]
  pseudobulk_qtnorm <- melt(pseudobulk_qtnorm, id.vars = c("sample_id", "id", "time"), variable.name = "gene")
  return(pseudobulk_qtnorm[order(sample_id),])
}

## center and scale within each gene, across all samples
standardize <- function(pseudobulk) {
  pseudobulk_standardized <- pseudobulk[, lapply(.SD, function(x){x-mean(x, na.rm = TRUE)/sd(x, na.rm = TRUE)}), by = gene, .SDcols = 'value']
  index <- pseudobulk[, .SD, by = gene, .SDcols = c('sample_id', 'id', 'time')]
  pseudobulk_standardized <- pseudobulk_standardized[, c('sample_id', 'id', 'time') := list(index$sample_id, index$id, index$time)]
  rm(index)
  pseudobulk_standardized <- pseudobulk_standardized[order(id, time, gene)]
  return(pseudobulk_standardized)
}

## get expression PCs
# `filtered` is just for naming plots
get_expr_pc <- function(pseudobulk, num_pc, filtered) {
  # wide format for pca, note the first column is sample_id, since data.table doesn't support rownames
  pseudobulk_pca_wf <- dcast(pseudobulk, sample_id + id + time ~ gene, value.var = "value")
  pca_index <- pseudobulk_pca_wf[, c('sample_id', 'id', 'time')]
  pseudobulk_pca_wf <- scale(pseudobulk_pca_wf[, !c("sample_id", "id", "time")])
  expr_usv <- svd(pseudobulk_pca_wf)
  pc_percent <- ((expr_usv$d^2)/sum((expr_usv$d^2)))*100
  # plot variance explained and pc1 vs pc2
  png(paste0(fig_dir, "exprPC_scree_lineage", which_data, "_", which_time, "_", filtered, ".png"))
  plot(pc_percent, pch = 20, xlim = c(1, num_pc), xlab = "PC Number", 
       ylab = "% Variance explained", main = "Scree plot")
  dev.off()
  pca_index[, c("PC1", "PC2") := list(expr_usv$u[, 1], expr_usv$u[, 2])]
  p <- ggplot(pca_index, aes(x = PC1, y = PC2, color = time)) + geom_point()
  ggsave(plot = p, filename = paste0(fig_dir, "pc1vs2_lineage", which_data, "_", which_time, "_", filtered, ".png"))
  expr_pc <- expr_usv$u[, 1:num_pc]
  rownames(expr_pc) <- pca_index$sample_id
  return(expr_pc)
}

## get cell line PCs
get_cellline_pc <- function(pseudobulk, num_pc, filtered) {
  # drop the genes that have zero variance in all the remaining cell lines for any time (if any), 
  # otherwise we can't standardize that column
  gene_day_sd <- pseudobulk %>% group_by(gene, time) %>% summarise(sd = sd(value))
  gene_day_sd <- gene_day_sd %>% group_by(gene) %>% filter(all(sd > 0))
  genes_clpca <- unique(gene_day_sd$gene)
  pseudobulk <- pseudobulk[gene %in% genes_clpca,]
  
  # wide format for cell line PCA
  pseudobulk_clpca_wf <- dcast(pseudobulk, id ~ time + gene, value.var = "value")
  clpca_index <- pseudobulk_clpca_wf[, "id"]
  pseudobulk_clpca_wf <- scale(pseudobulk_clpca_wf[, !"id"])
  cellline_usv <- svd(pseudobulk_clpca_wf)
  clpc_percent <- ((cellline_usv$d^2)/sum((cellline_usv$d^2)))*100
  png(paste0(fig_dir, "celllinePC_scree_lineage", which_data, "_", which_time, ".png"))
  plot(clpc_percent, pch = 20, xlim = c(1,num_pc), xlab = "cell line PC Number", 
       ylab = "% Variance explained", main = "Scree plot")
  dev.off()
  cellline_pc <- cellline_usv$u[, 1:num_pc]
  rownames(cellline_pc) <- clpca_index$id
  return(cellline_pc)
}

###############################
## perform preprocessing
###############################

raw_pseudobulk <- get_raw_pseudobulk()
subset_results <- subset_samples(raw_pseudobulk)
pseudobulk <- subset_results$pb
missing_samples <- subset_results$missing_samples
# find the cell lines with missing samples
# drop them in downstream analysis for `pseudobulk_filtered`
missing_celllines <- unique(sapply(missing_samples, function(x){strsplit(x, '_')}[[1]][1]))
pseudobulk_filtered <- pseudobulk[!id %in% missing_celllines,]

pseudobulk <- quantile_normalize(log_normalize(pseudobulk))
pseudobulk_filtered <- quantile_normalize(log_normalize(pseudobulk_filtered))

# for the full pseudobulk, can only get expression pcs
# keep it in case we have a model that doesn't use cell line PCs
# the pseudobulk we used to compute PC is log normalized and quantile normalized
# standardization is done inside the PCA function 
# `uf` is unfiltered
expr_pc_pb  <- get_expr_pc(pseudobulk, num_expr_pc, "uf")
saveRDS(expr_pc_pb, paste0(data_dir, "exprPC", num_expr_pc, "_lineage", which_data, "_", which_time, "_uf.rds"))
# save the standardized data after computing 
pseudobulk_preprocessed <- standardize(pseudobulk)
saveRDS(pseudobulk_preprocessed, paste0(data_dir, "pseudobulk_preprocessed_lineage", which_data, "_", which_time, "_uf.rds"))

# for the cell line filtered pseudobulk, also compute cell line pcs
expr_pc_pbf <- get_expr_pc(pseudobulk_filtered, num_expr_pc, "f")
saveRDS(expr_pc_pbf, paste0(data_dir, "exprPC", num_expr_pc, "_lineage", which_data, "_", which_time, "_f.rds"))
cellline_pc <- get_cellline_pc(pseudobulk_filtered, num_cellline_pc)
saveRDS(cellline_pc, paste0(data_dir, "celllinePC", num_cellline_pc, "_lineage", which_data, "_", which_time, ".rds"))
pseudobulk_filtered_preprocessed <- standardize(pseudobulk_filtered)
saveRDS(pseudobulk_filtered_preprocessed, paste0(data_dir, "pseudobulk_preprocessed_lineage", which_data, "_", which_time, "_f.rds"))

