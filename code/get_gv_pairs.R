gene.locs <- readRDS("../data/biomart_genepos.rds")
snp.locs <- read.table("../data/snp_locations_05cut.txt")
gene.locs$chromosome_name <- factor(paste0("chr", gene.locs$chromosome_name), levels=levels(snp.locs$chr))

for (c in levels(gene.locs$chromosome_name)) {
  chr.locs <- subset(gene.locs, chromosome_name==c)
  for (g in chr.locs$hgnc_symbol) {
    tss <- chr.locs[chr.locs$hgnc_symbol==g,"transcription_start_site"]
    if (g == chr.locs$hgnc_symbol[1]) {
      gv <- subset(snp.locs, (chr==c) & (pos >= tss-1e6) & (pos <= tss+1e6))
    } else {
      gv <- rbind(gv, subset(snp.locs, (chr==c) & (pos >= tss-1e6) & (pos <= tss+1e6)))
    }
  }
  saveRDS(gv, paste0("../data/gene_snps/", c, ".rds"))
  rm(chr.locs, tss, gv)
}

