output_dir <- "results/dynamic_pb/"
which_data <- "all"
which_time <- "day"
chr <- 13

results <- readRDS(paste0(output_dir, "lineage", which_data, "_", which_time, "_chr", chr, ".rds"))

# bonferroni correction for all gene-variant pairs 
# TODO: should we do this for all chromosomes together?
results$bonf_pvalue <- p.adjust(results$`p-value`, method = "bonferroni")

# get top p-value for each gene
top_hits <- results[, lapply(.SD, min, na.rm = TRUE), by = gene]

# FDR on gene level
top_hits$bh_pvalue <- p.adjust(top_hits$bonf_pvalue, method = "BH")
