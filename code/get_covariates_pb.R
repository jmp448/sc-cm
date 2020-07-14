library(Matrix)

cell_types <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15", 
          "iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC")
for (d in cell_types) {
  data <- readRDS(paste0("../data/static_pb/", d, "/preprocessed.rds"))
  svdecomp <- svd(data, nv=10)
  
  # save the scree plot
  evals <- sapply(svdecomp[[1]], function(x){(x^2)})
  pve <- evals/sum(evals)*100
  png(paste0("../figs/pb/scree_", d, ".png"))
  plot(seq(length(pve)), pve, xlab="nPCs", ylab="PVE", main=d)
  dev.off()
  
  # save the pcs
  v <- svdecomp$v
  rownames(v) <- colnames(data)
  colnames(v) <- paste0("PC", seq(1:ncol(v)))
  saveRDS(t(v), paste0("../data/static_pb/", d, "/exp_pcs.rds"))
  for (npc in seq(1,9)) {
    write.table(t(v[,1:npc]), file=paste0("../data/static_pb/", d, "/", npc, "pcs.txt"),
                quote=F, sep="\t", row.names=T, col.names=T)
  }
}

# look at factor matrix
#u <- svdecomp[[2]]
#u1 <- u[,1]
#names(u1) <- rownames(data)
#names(sort(abs(u1), decreasing=T))[1:10]