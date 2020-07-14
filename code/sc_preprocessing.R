library(Seurat)
library(reshape)
library(Matrix)
library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

sc <- readRDS("../../sc3/data/mitofilt_singlets.umap.cc.rds")

types <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15",
           "iPSC", "mesoderm", "EMT", "progenitor", "cardiomes", "CM", "EPDC")
for (ct in types) {
  # copy subsetted genotype files over
  stopifnot(file.copy(paste0("../data/static_pb/", ct, "/genotypes_05cut.txt"),
            paste0("../data/static_sc/", ct, "/genotypes_05cut.txt"), overwrite=T))
  stopifnot(file.copy(paste0("../data/static_pb/", ct, "/snp_locations_05cut.txt"),
            paste0("../data/static_sc/", ct, "/snp_locations_05cut.txt"), overwrite=T))
  
  # copy subsetted gene files over
  stopifnot(file.copy(paste0("../data/static_pb/", ct, "/gene_locations.txt"),
                      paste0("../data/static_sc/", ct, "/gene_locations.txt"), overwrite=T))
  
  # get list of genes to be included
  genes = read_tsv(paste0("../data/static_sc/", ct, "/gene_locations.txt"))$gene
  # get list of variants to be included
  snps = read_tsv(paste0("../data/static_sc/", ct, "/genotypes_05cut.txt"), skip=1, col_names=F)$X1
  # get list of individuals to be included
  inds.here = colnames(read.delim(paste0("../data/static_sc/", ct, "/genotypes_05cut.txt"), nrows=1))
  
  # subset inds and counts matrices
  if (substr(ct, 1, 3)=="day") {
    sc.ct = subset(sc, (demux.sng1 %in% inds.here) & (diffday==ct))
  } else {
    sc.ct = subset(sc, (demux.sng1 %in% inds.here) & (cell.type==ct))
  }
  sc.ct = sc.ct[genes,]
  
  # normalize, center, and scale
  sc.ct = NormalizeData(sc.ct)
  sc.ct = ScaleData(sc.ct, features=rownames(sc.ct))
  
  # save expression data
  counts.ct = sc.ct@assays[['RNA']]@scale.data
  saveRDS(counts.ct, file=paste0("../data/static_sc/", ct, "/exp.rds"))
  
  # save individual data
  inds.ct = as_tibble(sc.ct$demux.sng1)
  colnames(inds.ct) = "ind"
  write_tsv(inds.ct, paste0("../data/static_sc/", ct, "/inds.tsv"))
  
  # run PCA and save pcs
  sc.ct = RunPCA(sc.ct, npcs=100)
  pcs = sc.ct@reductions[['pca']]@cell.embeddings
  saveRDS(pcs, file=paste0("../data/static_sc/", ct, "/pcs.rds"))
  rm(sc.ct, inds.ct, pcs)
  
  # subset tests file to only those tests that will be performed for this cell type
  for (c in 1:22) {
    tests = read_tsv(paste0("../data/gene_snps/chr", c, ".tsv"))
    tests = filter(tests, (gene %in% genes) & (snp %in% snps))
    write_tsv(tests, paste0("../data/static_sc/", ct, "/gene_snps/chr", c, ".tsv"))
  }
}