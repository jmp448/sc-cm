library(lme4)
library(pbkrtest)

data_dir <- "data/dynamic_pb/"
output_dir <- "results/dynamic_pb/"
which_data <- "all"
which_time <- "day"
chr <- 13
num_expr_pc <- 5
num_cellline_pc <- 5

data <- readRDS(paste0(data_dir, "input_lineage", which_data, "_", which_time, "_chr", chr, ".rds"))

# create a vector of names for the PCs
expr_pcs <- c()
if (num_expr_pc > 0) {for (i in 1:num_expr_pc) {expr_pcs[i] <- paste0("exprPC", i)}}
cellline_pcs <- c()
if (num_cellline_pc > 0) {for (i in 1:num_cellline_pc) {cellline_pcs[i] <- paste0("celllinePC", i)}}

test_results <- list()
index <- 1
for (tid in unique(data$test_id)) {
  variable_columns <- c('value', 'time', 'genotype', 'id')
  variable_columns <- append(variable_columns, expr_pcs)
  variable_columns <- append(variable_columns, cellline_pcs)
  test_input <- data[test_id == tid, ..variable_columns]
  test_input$genotype_x_time <- test_input$genotype * test_input$time
  # compute the interaction between time and PCs
  for (i in 1:length(expr_pcs)) {
    colname <- paste0("exprPC", i, "_x_", "time")
    new_col <- test_input[, get(paste0("exprPC", i))] * test_input$time
    test_input[, eval(colname) := new_col]
  }
  for (i in 1:length(cellline_pcs)) {
    colname <- paste0("celllinePC", i, "_x_", "time")
    new_col <- test_input[, get(paste0("celllinePC", i))] * test_input$time
    test_input[, eval(colname) := new_col]
  }
  model <- lmer(value ~ (. - id) + (1 | id), data = test_input)
  coefs <- data.frame(coef(summary(model)))
  # use the Kenward-Roger approximation to get approximate degrees of freedom and the t-distribution to get p-values
  # source: https://www.r-bloggers.com/three-ways-to-get-parameter-specific-p-values-from-lmer/
  df.KR <- get_Lb_ddf(model, fixef(model))
  coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))  
  test_results[[index]] <- coefs
  print(paste0(index, " in ", length(unique(data$test_id))))
  index <- index + 1
}
names(test_results) <- unique(data$test_id)
saveRDS(test_results, paste0(output_dir, "lineage", which_data, "_", which_time, "_chr", chr, ".rds"))

