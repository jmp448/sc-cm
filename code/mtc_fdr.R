library(readr)
library(tibble)

cell_types <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15",
                "iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC")

for (d in cell_types) {
  for (npc in seq(0,9)) {
    cisqtl <- read_tsv(paste0("../data/", d, "/", npc, "_cisqtl.txt"))
    gene_ct <- table(cisqtl$gene)
    
    # bonferroni MTC
    cisqtl$bonf.p <- sapply(gene_ct[cisqtl$gene] * cisqtl$`p-value`, function(x){min(1,x)})
    
    # get the top pvalue for each gene
    cisqtl <- cisqtl[order(cisqtl$bonf.p),]
    cisqtl <- cisqtl[!duplicated(cisqtl$gene),]
    
    # FDR
    bh <- seq(nrow(cisqtl))/nrow(cisqtl)*0.05
    cisqtl <- cisqtl[1:max(which(cisqtl$bonf.p<bh)),]
    saveRDS(cisqtl, paste0("../data/", d, "/", npc, "_sighits.txt"))
  }
}
