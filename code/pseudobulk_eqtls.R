library(MatrixEQTL)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

cell_types <- c("day0", "day1", "day3", "day5", "day7", "day11", "day15",
                "iPSC", "mesoderm", "EMT", "cardiomes", "progenitor", "CM", "EPDC")

for (d in cell_types) {
  for (npc in c(0, 1,2,3,4,5,6,7,8,9)) {
    # Genotype file name
    SNP_file_name = paste0("../data/static_pb/", d, "/genotypes_05cut.txt")
    snps_location_file_name = paste0("../data/static_pb/", d, "/snp_locations_05cut.txt")
    
    # Gene expression file name
    expression_file_name = paste0("../data/static_pb/", d, "/pseudobulk_preprocessed.txt")
    gene_location_file_name = paste0("../data/static_pb/", d, "/gene_locations.txt")
    
    # Covariates file name
    if (npc==0) {
      covariates_file_name = character()
    } else {
      covariates_file_name = paste0("../data/static_pb/", d, "/", npc, "pcs.txt")
    }
    
    # Output file name
    output_file_name_cis = paste0("../results/static_pb/", d, "/", npc, "_cisqtl.txt")
    output_file_name_tra = tempfile()
    
    # Only associations significant at this level will be saved
    pvOutputThreshold_cis = 1;
    pvOutputThreshold_tra = 0;
    
    # Error covariance matrix
    # Set to numeric() for identity.
    errorCovariance = numeric();
    
    # Distance for local gene-SNP pairs
    cisDist = 1e6;
    
    ## Load genotype data
    
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name);
    
    ## Load gene expression data
    
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";      # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## Load covariates
    
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
      cvrt$LoadFile(covariates_file_name);
    }
    
    ## Run the analysis
    snpspos = read.table(snps_location_file_name, header = TRUE, row.names=NULL, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name     = output_file_name_tra,
      pvOutputThreshold     = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = TRUE);
    
    unlink(output_file_name_tra);
    
    ## Results:
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  }
}
