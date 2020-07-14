library(readr)
library(dplyr)
library(tibble)

gene.locs <- tibble(readRDS("../data/biomart_genepos.rds"))
snp.locs <- read_tsv("../data/snp_locations_05cut.txt", skip=1, col_names=F)
colnames(snp.locs) <- c("snp", "chr", "pos")
gene.locs$chromosome_name <- paste0("chr", gene.locs$chromosome_name)

get_snps <- function(g, chrom, slocs=snp.locs, glocs=gene.locs) {
  tss = filter(glocs, hgnc_symbol==g)$transcription_start_site
  gv = select(filter(slocs, (chr==chrom) & (pos >= tss-1e6) & (pos <= tss+1e6)), snp)
  gv$gene = g
  gv
}

for (c in unique(gene.locs$chromosome_name)) {
  chr.locs <- filter(gene.locs, chromosome_name==c)
  tests.c <- bind_rows(lapply(chr.locs$hgnc_symbol, get_snps, chrom=c))
  write_tsv(tests.c, paste0("../data/gene_snps/", c, ".tsv"))
}
