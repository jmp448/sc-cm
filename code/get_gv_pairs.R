gene.locs <- readRDS("data/biomart_genepos.rds")
snp.locs <- read.table("data/snp_locations_05cut.txt")
gene.locs$chromosome_name <- factor(paste0("chr", gene.locs$chromosome_name), levels=levels(snp.locs$chr))

for (c in levels(gene.locs$chromosome_name)) {
  chr.locs <- subset(gene.locs, chromosome_name==c)
  for (g in chr.locs$hgnc_symbol) {
    tss <- chr.locs[chr.locs$hgnc_symbol==g,"transcription_start_site"]
    gv <- subset(snp.locs, (chr==c) & (pos >= tss-1e6) & (pos <= tss+1e6))
    if (dim(gv)[1] == 0) {next}
    # the gene associated with this variant
    gv$gene <- g
    # keep the variant name for easier access of genotype data
    gv$variant <- rownames(gv)
    if (g == chr.locs$hgnc_symbol[1]) {
      gv_all <- gv
    } else {
      gv_all <- rbind(gv_all, gv)
    }
  }
  # since if a variant appears for multiple genes, the rowname wouldn't make sense
  # reset it to be clear
  rownames(gv_all) <- NULL
  saveRDS(gv_all, paste0("data/gene_snps/", c, ".rds"))
  rm(chr.locs, tss, gv, gv_all)
}

