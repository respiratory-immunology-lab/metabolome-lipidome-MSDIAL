add_hmdb <- function(metab_SE, hmdb, mass_tol, cores) {
  pkgs <- c('foreach', 'doSNOW', 'itertools', 'dplyr', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # Set ion mass
  ion <- 1.007276
  # Transform everything into a vector for faster looping
  hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
  hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  Mass_db <- as.vector(as.numeric(hmdb$monisotopic_molecular_weight))
  KEGG_db <- as.vector(hmdb$kegg)
  Name_db <- as.vector(hmdb$name)
  HMDB_id_db <- as.vector(hmdb$accession)
  Mass_data <- as.vector(rowData(metab_SE)$`info.Average Mz`)
  Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  HMDB_data <- rep(NA, length(Mass_data))
  HMDB_id <- rep(NA, length(Mass_data))
  KEGG_data <- rep(NA, length(Mass_data))
  # Get masses corrected for ion precursors
  for (m in 1:length(Mass_data)){
    if (Adduct_data[m] %in% c('[M+2H]2+','[M+H]+','[M+Na]+','[M+NH4]+','[M+NH4]2+')){
      Mass_data[m] <- Mass_data[m]-ion}
    else {
      Mass_data[m] <- Mass_data[m]+ion}}
  
  Matrix <- matrix(nrow = length(Mass_db), ncol = length(Mass_data))
  rownames(Matrix) <- Mass_db
  colnames(Matrix) <- Mass_data
  #Setup clusters
  cores=cores
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  message("HMDB Annotation Starting - Get a Coffee While You Wait :)")
  start_time <- Sys.time()
  TF_DF <- foreach(i = isplitCols(Matrix, chunks=cores), .combine = "cbind", .packages = c("dplyr")#, .options.snow=opts
  ) %dopar% {
    for (n in 1:nrow(i)) {
      for (m in 1:ncol(i)) {
        i[n,m] <- between(as.numeric(rownames(i)[n]), as.numeric(colnames(i)[m])-mass_tol, as.numeric(colnames(i)[m])+mass_tol)
      }
    }
    i
  }
  end_time <- Sys.time()
  message("Minutes Taken: ", round(end_time - start_time,2))
  
  stopCluster(cl)
  
  #Make names
  for (i in c(1:ncol(TF_DF))) {
    if (sum(TF_DF[,i]*1) > 0) {
      index = which(TF_DF[,i] == T) %>% as.numeric()
      HMDB_data[i] <- paste0(Name_db[index], collapse = ';')
      HMDB_id[i] <- paste0(Name_db[index], collapse = ';')
      KEGG_data[i] <- paste0(Name_db[index], collapse = ';')
    }
  }
  # Add new information to SE experiment object
  rowData(metab_SE)$HMDB <- HMDB_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$KEGG <- KEGG_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$HMDB_accession <- HMDB_id %>% dplyr::na_if(.,'')
  return(metab_SE)
}