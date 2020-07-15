data_dir <- "data/dynamic_pb/"
chr <- "13"
num_expr_pc <- 10
num_cellline_pc <- 10
which_data <- "all"
which_time <- "day"

expr <- readRDS(paste0(data_dir, "pseudobulk_preprocessed_lineage", which_data, "_", which_time, "_", "f.rds"))
expr_pc <- readRDS(paste0(data_dir, "exprPC30_lineage", which_data, "_", which_time, "_", "f.rds"))
cellline_pc <- readRDS(paste0(data_dir, "celllinePC10_lineage", which_data, "_", which_time, ".rds"))
genotypes <- read.table("data/genotypes_05cut.txt", header=T, row.names=1)

expr_pc <- expr_pc[, 1:num_expr_pc]
expr_pc_names <- c()
for (i in 1:num_expr_pc) {expr_pc_names[i] <- paste0("expr_PC", i)}
colnames(expr_pc) <- expr_pc_names
expr_pc <- data.table(sample_id = rownames(expr_pc), expr_pc)

cellline_pc <- cellline_pc[, 1:num_cellline_pc]
cellline_pc_names <- c()
for (i in 1:num_cellline_pc) {cellline_pc_names[i] <- paste0("cellline_PC", i)}
colnames(cellline_pc) <- cellline_pc_names
cellline_pc <- data.table(id = rownames(cellline_pc), cellline_pc)

gene_snps <- readRDS(paste0("data/gene_snps/chr", chr, ".rds"))
# each line of the output is a test to run, containing all the information for that test
# looping through all gene-variant pair
for (i in 1:dim(gene_snps)[1]) {
  g <- gene_snps[i, "gene"]
  v <- gene_snps[i, "variant"]
  g_expr <- expr[gene == g]
  g_expr$variant <- v
  
  # use variant and id to get genotype
  v_geno <- t(genotypes[v,])
  v_geno <- data.table(id = rownames(v_geno), genotype = c(v_geno))
  g_expr <- merge(g_expr, v_geno, by = "id")
  
  # include expression PCs and cell line PCs
  g_expr <- merge(g_expr, expr_pc, by = "sample_id")
  g_expr <- merge(g_expr, cellline_pc, by = "id")
  
  if (i == 1) {
    output <- g_expr
  } else {
    output <- rbind(output, g_expr)
  }
  print(paste0(i, " in ", dim(gene_snps)[1]))
}
saveRDS(output, paste0(data_dir, "input_lineage", which_data, "_", which_time, ".rds"))



