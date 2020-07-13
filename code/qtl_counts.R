library(ggplot2)
library(readr)
library(reshape2)
library(corrplot)
library(stats)

cell_types <- c("iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC")

for (ct in cell_types) {
  for (npc in seq(0,9)) {
    sighits <- readRDS(paste0("../data/", ct, "/", npc, "_sighits.txt"))
    if ((ct=="iPSC") & (npc==0)) {
      counts <- data.frame("type"="iPSC", "pcs"=0, "egenes"=nrow(sighits))
    } else {
      counts <- rbind(counts, data.frame("type"=ct, "pcs"=npc, "egenes"=nrow(sighits)))
    }
  }
}
counts$pcs <- factor(counts$pcs)
counts$type <- factor(counts$type)

ggplot(data=counts, aes(x=type, y=egenes, fill=pcs)) +
  geom_bar(stat="identity", position=position_dodge())

days <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15")
for (d in days) {
  for (npc in seq(0,9)) {
    sighits <- readRDS(paste0("../data/", d, "/", npc, "_sighits.txt"))
    if ((d=="day0") & (npc==0)) {
      counts <- data.frame("day"="day0", "pcs"=0, "egenes"=nrow(sighits))
    } else {
      counts <- rbind(counts, data.frame("day"=d, "pcs"=npc, "egenes"=nrow(sighits)))
    }
  }
}
counts$pcs <- factor(counts$pcs)
counts$day <- factor(counts$day)

ggplot(data=counts, aes(x=day, y=egenes, fill=pcs)) +
  geom_bar(stat="identity", position=position_dodge())

# compare the pvals
cell_types <- c("iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC",
                "day0", "day1", "day3", "day5", "day7", "day11", "day15")
sighits <- readRDS("../data/iPSC/0_sighits.txt")
sighits$"type" <- "iPSC"
for (ct in cell_types[-c(1)]) {
  sighits_i <- readRDS(paste0("../data/", ct, "/0_sighits.txt"))
  sighits_i$"type" <- ct
  sighits <- rbind(sighits, sighits_i)
  rm(sighits_i)
}

# get list of unique egene - evariant pairs
sighits$gv <- mapply(paste, sighits$gene, sighits$SNP, MoreArgs=list(sep="_"))
# egenes <- unique(sighits$gene)
# evariants <- unique(sighits$SNP)
sig_gvs <- unique(sighits$gv)

# get the full eQTL analysis results
cisqtl <- read_tsv(paste0("../data/iPSC/0_cisqtl.txt"))
gene_ct <- table(cisqtl$gene)
cisqtl$gv <- mapply(paste, cisqtl$gene, cisqtl$SNP, MoreArgs=list(sep="_"))
cisqtl <- subset(cisqtl, gv %in% sig_gvs)
# cisqtl <- subset(cisqtl, gene %in% egenes)
cisqtl$bonf.p <- sapply(gene_ct[cisqtl$gene] * cisqtl$`p-value`, function(x){min(1,x)})
cisqtl <- subset(cisqtl, select=c("SNP", "gene", "p-value", "bonf.p"))
cisqtl$type <- "iPSC"

for (d in cell_types[-c(1)]) {
  cisqtl_i <- read_tsv(paste0("../data/", d, "/0_cisqtl.txt"))
  gene_ct <- table(cisqtl$gene)
  cisqtl_i$gv <- mapply(paste, cisqtl_i$gene, cisqtl_i$SNP, MoreArgs=list(sep="_"))
  cisqtl_i <- subset(cisqtl_i, gv %in% sig_gvs)
  # cisqtl_i <- subset(cisqtl_i, gene %in% egenes)
  cisqtl_i$bonf.p <- sapply(gene_ct[cisqtl_i$gene] * cisqtl_i$`p-value`, function(x){min(1,x)})
  cisqtl_i <- subset(cisqtl_i, select=c("SNP", "gene", "p-value", "bonf.p"))
  cisqtl_i$type <- d
  cisqtl <- rbind(cisqtl, cisqtl_i)
  rm(cisqtl_i)
}

# subset to those SNP-gene pairs that appear in all cell
# cisqtl$gv <- mapply(paste, cisqtl$gene, cisqtl$SNP, MoreArgs=list(sep="_"))
# cisqtl <- subset(cisqtl, gv %in% names(table(cisqtl$gv)[table(cisqtl$gv)==14]))
# cisqtl <- subset(cisqtl, select=c("SNP", "gene", "bonf.p", "type"))
# cisqtl <- cisqtl[order(cisqtl$type, cisqtl$gene, cisqtl$SNP),]
# duos <- nrow(cisqtl)/14
# stopifnot(all.equal(cisqtl$SNP[1:duos], cisqtl$SNP[(duos+1):(2*duos)]))

# reshape so each type gets its own column
celltypes <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15",
               "iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC")
cisqtl$type <- factor(cisqtl$type, levels=celltypes)
pvals <- dcast(cisqtl, SNP + gene ~ type, value.var="bonf.p")
bonf <- as.matrix(pvals[,3:16])
bonf <- -log10(bonf)
cor.p <- cor(bonf, use="complete.obs")
corrplot.mixed(cor.p, tl.pos="lt", number.cex=0.7, lower.col="black")

# how does correlation of expression look?
for (ct in cell_types) {
  if (ct == "iPSC") {
    exp <- read.table(paste0("../data/", ct, "/pseudobulk_lognorm.txt"))
    exp <- data.frame(apply(exp, 1, mean))
  } else {
    exp_ct <- read.table(paste0("../data/", ct, "/pseudobulk_lognorm.txt"))
    genes <- intersect(rownames(exp), rownames(exp_ct))
    exp <- exp[genes,]
    exp_ct <- exp_ct[genes,]
    exp_ct <- data.frame(apply(exp_ct, 1, mean))
    exp <- cbind(exp, exp_ct)
  }
}
colnames(exp) <- cell_types
# pseudobulk expression correlation
corrplot(cor(exp), method="number")

# is it too good to be true? look at marker genes
annotation <- data.frame(
  x = c(exp['TNNT2', 'iPSC'], exp['POU5F1', 'iPSC']),
  y = c(exp['TNNT2', 'CM'], exp['POU5F1', 'CM']),
  label = c("TNNT2", "POU5F1")
)

ggplot(exp, aes(x=iPSC, y=CM)) +
  geom_point() +
  geom_label(data=annotation, aes(x=x, y=y, label=label))

# how about the scatter plots among -log10 bonf pvals
log10.pvals <- cisqtl
log10.pvals$bonf.p <- -log10(log10.pvals$bonf.p)
log10.pvals <- dcast(log10.pvals, SNP + gene ~ type, value.var="bonf.p")
ggplot(log10.pvals, aes(x=iPSC, y=CM)) +
  geom_point()