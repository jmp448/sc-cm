library(Seurat)
library(reshape)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(reshape2)

sc <- readRDS("../data/mitofilt_singlets.umap.cc.rds")

# Tweak clustering
sc <- FindClusters(sc, resolution=0.2)
DimPlot(sc, group.by="seurat_clusters")
saveRDS(sc, "../data/mitofilt_singlets.umap.cc.rds")

# Cell typing
FeaturePlot(sc, features="TNNT2")  # cluster 5 CMs
FeaturePlot(sc, features="DCN") # cluster 1 EPDCs
FeaturePlot(sc, features="POU5F1") # cluster 2 iPSCs
FeaturePlot(sc, features="TBXT") #cluster 0 mesoderm
FeaturePlot(sc, features="MESP1") #cluster 4 cardiac mesoderm
FeaturePlot(sc, features="KDR") # cluster 3 progenitor cells?
FeaturePlot(sc, features="VIM") # cluster 6 EMT

celltype <- function(clust) {
  if (clust==0) {
    return("mesoderm")
  } else if (clust==1) {
    return("EPDC")
  } else if (clust==2) {
    return("iPSC")
  } else if (clust==3) {
    return("progenitor")
  } else if (clust==4) {
    return("cardiomes")
  } else if (clust==5) {
    return("CM")
  } else if (clust==6) {
    return("EMT")
  }
}

# clean up the seurat object
ctlevs <- c("iPSC", "mesoderm", "EMT", "cardiomes",
            "progenitor", "CM", "EPDC")
sc$cell.type <- factor(sapply(sc$seurat_clusters, celltype), levels=ctlevs)

clean_diffday <- function(dday) {
  dday <- gsub("Day ", "day", dday)
  dday
}
daylevs <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15")
sc$diffday <- factor(sapply(sc$diffday, clean_diffday), levels=daylevs)
sc$demux.sng1 <- as.factor(sc$demux.sng1)

# Get important data
counts <- sc@assays[["RNA"]]@counts
type <- sc$cell.type
day <- sc$diffday
inds <- sc$demux.sng1

# subset counts to what appears in biomart
# this section can be removed after full seurat obj use
hgnc <- readRDS("../../cm_eqtl/data/biomart_genepos.rds")
hgnc <- subset(hgnc, hgnc_symbol %in% rownames(counts))
counts <- counts[rownames(counts) %in% hgnc$hgnc_symbol,]

# get snp and genotype matrices
snps <- read.table("../../cm_eqtl/data/genotypes_05cut.txt", header=T, row.names=1)
snp_locs <- read.table("../../cm_eqtl/data/snp_locations_05cut.txt", header=T, row.names=1)

# get pseudobulk for each desired combo
getPseudobulk <- function(mat, ind) {
  mat.summary <- do.call(cbind, lapply(levels(ind), function(i) {
    cells <- names(ind)[ind==i]
    pseudobulk <- rowSums(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(ind)
  return(mat.summary)
}


#should we log transform the data? does log transformation make it look more gaussian?
#exp <- as.vector(mat.preprocessed)
#ggplot(data.frame(exp), aes(x=exp)) + geom_density()
# it seems to help

## look at the representation of individuals for each cell type/ day
df <- data.frame(inds, day, type)
# ggplot(df) + geom_bar(aes(x=type))
# ggplot(df) + geom_bar(aes(x=type, fill=inds))
ggplot(df) + geom_bar(aes(x=inds)) + facet_grid(type) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

# ggplot(df) + geom_bar(aes(x=day))
# ggplot(df) + geom_bar(aes(x=day, fill=inds))
ggplot(df) + geom_bar(aes(x=inds)) + facet_grid(day) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

# look at the distribution of read counts as well
for (ct in levels(type)) {
  mat <- counts[,names(type)[type==ct]]
  inds_ct <- inds[names(type)[type==ct]]
  
  # raw pseudobulk
  mat.summary <- getPseudobulk(mat, inds_ct)
  mat.summary <- colSums(mat.summary)
  if (ct == "iPSC") {
    df2 <- data.frame(t(mat.summary))
    df2$"type" <- "iPSC"
  } else {
    df3 <- data.frame(t(mat.summary))
    df3$"type" <- ct
    df2 <- rbind(df2, df3)
  }
}
df2 <- melt(df2, id.vars=c("type"), variable.name="individual",
            value.name="readcounts")
df2$type <- factor(df2$type)
ggplot(df2, aes(x=individual, y=readcounts)) + geom_bar(stat="identity") + 
  facet_grid(rows=df2$type) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

for (d in levels(day)) {
  mat <- counts[,names(day)[day==d]]
  inds_d <- inds[names(day)[day==d]]
  
  # raw pseudobulk
  mat.summary <- getPseudobulk(mat, inds_d)
  mat.summary <- colSums(mat.summary)
  if (d == "day0") {
    df2 <- data.frame(t(mat.summary))
    df2$"day" <- "day0"
  } else {
    df3 <- data.frame(t(mat.summary))
    df3$"day" <- d
    df2 <- rbind(df2, df3)
  }
}
df2 <- melt(df2, id.vars=c("day"), variable.name="individual",
            value.name="readcounts")
df2$day <- factor(df2$day, levels=paste0("day", c(0,1,3,5,7,11,15)))
ggplot(df2, aes(x=individual, y=readcounts)) + geom_bar(stat="identity") + 
  facet_grid(rows=df2$day) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
