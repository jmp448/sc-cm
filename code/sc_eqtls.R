library(nlme)
library(pbkrtest)

inputArgs <-  commandArgs(TRUE)
cell_type <- inputArgs[1]
chrom <- inputArgs[2]
npcs <- inputArgs[3]

expr <- readRDS(paste0("../data/sc/", cell_type, "/expr.rds"))
tests <- readRDS(paste0("../data/gene_snps/chr", chrom, ".rds"))
genes <- tests$gene
genotypes <- readRDS(paste0("../data/sc/", cell_type, "/"))
inds <- readRDS()


eQTLtest <- function(gene, SNP) {
  exp <- data.frame(expr[,g])
  inds <- data.frame()
}
for(i in c(start:end)){
  print(i)
  geno_subset <- geno[ ,snps.select[[i]][ ,3]]
  geno_subset <- as.data.frame(geno_subset)
  print(dim(geno_subset))
  if(dim(geno_subset)[2] > 0){
    for(igeno in c(1:dim(geno_subset)[2])){
      #print(igeno)
      fixed_eff <- cbind(geno_subset[ ,igeno], expr_PCs[ ,1:npcs], geno_PCs[ ,1:5])
      colnames(fixed_eff)[1] <- colnames(geno_subset)[igeno]
      data_df <- data.frame(expr[ ,i], fixed_eff, sample_assignments$individuals)
      colnames(data_df) <- c("expression", colnames(fixed_eff), "individuals")
      eQTL.model <- lmer(expression ~ (. - individuals) + (1 | individuals), data = data_df)
      icoefs <- data.frame(coef(summary(eQTL.model)))
      df.KR <- get_Lb_ddf(eQTL.model, fixef(eQTL.model))
      icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
      if(igeno == 1){
        coef_snps <- icoefs[2, ]
      }else{
        coef_snps <- rbind(coef_snps, icoefs[2, ])
      }
    }
    coef_snps_by_gene[[index]] <- coef_snps
  }else{
    coef_snps_by_gene[[index]] <- NA
  }
  index <- index + 1
}

names(coef_snps_by_gene) <- names(snps.select)[start:end]
saveRDS(coef_snps_by_gene, paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/coef_snps_by_gene_", start, "to", end, ".rds"))