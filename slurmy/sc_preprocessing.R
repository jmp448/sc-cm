library(Seurat)
library(reshape)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(reshape2)

sc <- readRDS("../../sc3/data/mitofilt_singlets.umap.cc.rds")

counts <- sc@assays[["RNA"]]@counts
type <- sc$cell.type
day <- sc$diffday
inds <- sc$demux.sng1

# subset counts to what appears in biomart
# this section can be removed after full seurat obj use
hgnc <- readRDS("../data/biomart_genepos.rds")
hgnc <- subset(hgnc, hgnc_symbol %in% rownames(counts))
counts <- counts[rownames(counts) %in% hgnc$hgnc_symbol,]

# get snp and genotype matrices
snps <- read.table("../../cm_eqtl/data/genotypes_05cut.txt", header=T, row.names=1)
snp_locs <- read.table("../../cm_eqtl/data/snp_locations_05cut.txt", header=T, row.names=1)

## call function for each cell type
for (ct in levels(type)) {
  mat <- counts[,names(type)[type==ct]]
  inds_ct <- inds[names(type)[type==ct]]
  
  saveRDS(mat, paste0("../data/static_sc/", ct, "/raw.rds"))
  
  # subset to individuals that have at least 100k UMI
  reads <- colSums(mat)
  inds.here <- names(reads[reads>100000])
  mat.summary <- mat.summary[,inds.here]
  snps_ct <- snps[,inds.here]
  write.table(snps_ct, paste0("../data/static_sc/", ct, "/genotypes.txt"),
              quote=F, sep="\t")
  
  # subset to snps with MAF>=0.05 for the individuals present
  afs <- rowSums(snps_ct, na.rm=TRUE)/(2*ncol(snps_ct))
  mafs <- sapply(afs, function(x){min(x,1-x)})
  mafs <- mafs[mafs>=0.05]
  snps_ct <- snps_ct[names(mafs),]
  write.table(snps_ct, file=paste0("../data/static_sc/", ct, "/genotypes_05cut.txt"),
              quote=F, sep="\t")
  snp_ct_locs <- snp_locs[names(mafs),]
  write.table(snp_ct_locs, file=paste0("../data/static_sc/", ct, "/snp_locations_05cut.txt"),
              quote=F, sep="\t")
  
  # subset to genes that were previously evaluated
  n55.genes <- rownames(mat.summary[rowSums(mat.summary>=5)>=5,])
  mat.summary <- mat.summary[n55.genes,]
  # subset to genes w nonzero variance
  sd_ct <- apply(mat.summary, 1, sd)
  mat.summary <- mat.summary[sd_ct!=0,]
  
  # log normalize
  mat.lognorm <- apply(mat.summary, 2, function(c){log(10^4*c/sum(c)+1)})
  saveRDS(mat.lognorm, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/pseudobulk_lognorm.rds"))
  write.table(mat.lognorm, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/pseudobulk_lognorm.txt"),
              quote=F, sep="\t", row.names=T, col.names=T)
  
  # center and scale
  mat.preprocessed <- t(apply(mat.lognorm, 1, function(g){(g-mean(g))/sd(g)}))
  saveRDS(mat.preprocessed, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/pseudobulk_preprocessed.rds"))
  write.table(mat.preprocessed, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/pseudobulk_preprocessed.txt"),
              quote=F, sep="\t", row.names=T, col.names=T)
  
  # subset gene list to match those included
  hgnc_ct <- subset(hgnc, hgnc_symbol %in% rownames(mat.preprocessed))
  hgnc_ct <- hgnc_ct[match(rownames(mat.preprocessed), hgnc_ct$hgnc_symbol),]
  colnames(hgnc_ct) <- c("gene", "chr", "tss")
  hgnc_ct$end <- hgnc_ct$tss + 1
  hgnc_ct$chr <- paste0("chr", hgnc_ct$chr)
  saveRDS(hgnc_ct, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/gene_locations.rds"))
  write.table(hgnc_ct, paste0("../../cm_eqtl/data/pseudobulk/", ct, "/gene_locations.txt"),
              quote=F, sep="\t", row.names=F)
}
rm(hgnc_ct, snp_ct_locs, snps_ct, ct, inds_ct, sd_ct)

## call function for each day
for (d in levels(day)) {
  mat <- counts[,names(day)[day==d]]
  inds_d <- inds[names(day)[day==d]]
  
  # raw pseudobulk
  mat.summary <- getPseudobulk(mat, inds_d)
  saveRDS(mat.summary, paste0("../../cm_eqtl/data/pseudobulk/", d, "/pseudobulk_raw.rds"))
  
  # subset to individuals that have at least 100k UMI
  reads <- colSums(mat.summary)
  inds.here <- names(reads[reads>100000])
  mat.summary <- mat.summary[,inds.here]
  snps_d <- snps[,inds.here]
  write.table(snps_d, paste0("../../cm_eqtl/data/pseudobulk/", d, "/genotypes.txt"),
              quote=F, sep="\t")
  
  # subset to snps with MAF>=0.05 for the individuals present
  afs <- rowSums(snps_d, na.rm=TRUE)/(2*ncol(snps_d))
  mafs <- sapply(afs, function(x){min(x,1-x)})
  mafs <- mafs[mafs>=0.05]
  snps_d <- snps_d[names(mafs),]
  write.table(snps_d, file=paste0("../../cm_eqtl/data/pseudobulk/", d, "/genotypes_05cut.txt"),
              quote=F, sep="\t")
  snp_d_locs <- snp_locs[names(mafs),]
  write.table(snp_d_locs, file=paste0("../../cm_eqtl/data/pseudobulk/", d, "/snp_locations_05cut.txt"),
              quote=F, sep="\t")
  
  # subset to genes with at least 5 reads in at least 5 samples
  n55.genes <- rownames(mat.summary[rowSums(mat.summary>=5)>=5,])
  mat.summary <- mat.summary[n55.genes,]
  # subset to genes w nonzero variance
  sd_d <- apply(mat.summary, 1, sd)
  mat.summary <- mat.summary[sd_d!=0,]
  
  # log normalize
  mat.lognorm <- apply(mat.summary, 2, function(c){log(10^4*c/sum(c)+1)})
  saveRDS(mat.lognorm, paste0("../../cm_eqtl/data/pseudobulk/", d, "/pseudobulk_lognorm.rds"))
  write.table(mat.lognorm, paste0("../../cm_eqtl/data/pseudobulk/", d, "/pseudobulk_lognorm.txt"),
              quote=F, sep="\t", row.names=T, col.names=T)
  
  # center and scale
  mat.preprocessed <- t(apply(mat.lognorm, 1, function(g){(g-mean(g))/sd(g)}))
  saveRDS(mat.preprocessed, paste0("../../cm_eqtl/data/pseudobulk/", d, "/pseudobulk_preprocessed.rds"))
  write.table(mat.preprocessed, paste0("../../cm_eqtl/data/pseudobulk/", d, "/pseudobulk_preprocessed.txt"),
              quote=F, sep="\t", row.names=T, col.names=T)
  
  # subset gene list to match those included
  hgnc_d <- subset(hgnc, hgnc_symbol %in% rownames(mat.preprocessed))
  hgnc_d <- hgnc_d[match(rownames(mat.preprocessed), hgnc_d$hgnc_symbol),]
  colnames(hgnc_d) <- c("gene", "chr", "tss")
  hgnc_d$end <- hgnc_d$tss + 1
  hgnc_d$chr <- paste0("chr", hgnc_d$chr)
  saveRDS(hgnc_d, paste0("../../cm_eqtl/data/pseudobulk/", d, "/gene_locations.rds"))
  write.table(hgnc_d, paste0("../../cm_eqtl/data/pseudobulk/", d, "/gene_locations.txt"),
              quote=F, sep="\t", row.names=F)
}