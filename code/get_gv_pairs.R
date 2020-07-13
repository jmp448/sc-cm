gene.locs <- readRDS("../data/biomart_genepos.rds")
gene.locs$chr <- factor(paste0("chr", gene.locs$chromosome_name))
rownames(gene.locs) <- gene.locs$hgnc_symbol
genes <- rownames(gene.locs)

snp.locs <- read.table("../data/snp_locations_05cut.txt")

# write a function to take each gene, find its location, find all snps w/in 1MB of the TSS
find.snps <- function(gene) {
  g.snps <- subset(snp.locs, (chr==gene.locs[gene, "chr"]) & (pos >= gene.locs[gene, "transcription_start_site"]-1e6) & (pos <= gene.locs[gene, "transcription_start_site"]+1e6))
  g.snps
}
g <- genes[1]
gv <- find.snps(g)
for (g in genes[2:length(genes)]) {
  gv <- rbind(gv, find.snps(g))
}