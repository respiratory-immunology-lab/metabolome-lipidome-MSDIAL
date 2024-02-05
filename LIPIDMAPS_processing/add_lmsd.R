add_lmsd <- function(metab_SE, lmsd, mass_tol = 0.002, cores = NA) {
  # Load required packages
  pkgs <- c('foreach', 'doParallel', 'data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'doSNOW')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # Set ion mass
  ion <- 1.007276
  lmsd$EXACT_MASS <- as.numeric(lmsd$EXACT_MASS)
  Mass_db <- lmsd$EXACT_MASS
  Mass_data <- data.frame(msdial_mz = rowData(metab_SE)$`info.Average Mz`)
  Mass_data$msdial_mz <- as.numeric(Mass_data$msdial_mz)
  rownames(Mass_data) <- rownames(rowData(metab_SE))
  Mass_data$Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  Mass_data$corrected_mz <- NA
  # Get masses corrected for ion precursors
  for (m in 1:nrow(Mass_data)){
    if (Mass_data$Adduct_data[m] %in% c("[M+H]+","[M+NH4]+","[M+2H]2+","[M+H-H2O]+", "[M+Na]+")){
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]-ion}
    else {
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]+ion
    }
  }
  #setup parallel backend to use many processors
  if(is.na(cores) | cores > detectCores()){
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
  } else {
    cl <- makeCluster(cores)
  }
  registerDoSNOW(cl)
  #Progress Bar
  iterations <- nrow(Mass_data)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  ## Compare masses corrected with lmsd "exact masses" to find annotations that are with tolerance range
  full_lmsd_ann <- foreach(n=1:nrow(Mass_data), .combine=rbind, .packages = "dplyr", .options.snow = opts) %dopar% {
    row_num <- c()
    for (m in 1:length(Mass_db)) {
      if(between(Mass_db[m], Mass_data$corrected_mz[n]-mass_tol, Mass_data$corrected_mz[n]+mass_tol)==TRUE)
        row_num <- c(row_num, m)
    }
    if(is.null(row_num)) { ## if no matches are found we want to add the data as NA
      temp <- data.frame(matrix(ncol = ncol(lmsd)))
      colnames(temp) <- colnames(lmsd)
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- 1 ## delta 1 = no match
    } else {
      temp <- lmsd[row_num,]
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- abs(Mass_data$corrected_mz[n] - lmsd$EXACT_MASS[row_num])
    }
    temp
  }
  close(pb)
  stopCluster(cl)
  
  # arrange smallest delta for each lipidID
  full_lmsd_ann <- full_lmsd_ann %>%
    arrange(delta, by_group = T) %>%
    arrange(factor(LipidID, levels = rownames(Mass_data))) # maintain original order
  
  # Remove duplicate matches for lipidID
  distinct_lmsd_ann <- full_lmsd_ann %>%
    distinct(LipidID, NAME, SYSTEMATIC_NAME, .keep_all=T)
  
  
  # Select columns of interest
  lmsd_ann_sub <- distinct_lmsd_ann[, c("LipidID","corrected_mz","delta",
                                        "NAME","SYSTEMATIC_NAME", "CATEGORY",
                                        "MAIN_CLASS","EXACT_MASS", "ABBREVIATION",
                                        "SYNONYMS","KEGG_ID","HMDB_ID", "SUB_CLASS",
                                        "CLASS_LEVEL4")]
  # Change column names
  colnames(lmsd_ann_sub) <- c("LipidID","corrected_mz","delta",
                              "LMSD_NAME","LMSD_SYSTEMATIC_NAME", "LMSD_CATEGORY",
                              "LMSD_MAIN_CLASS","LMSD_EXACT_MASS", "LMSD_ABBREVIATION",
                              "LMSD_SYNONYMS","LMSD_KEGG_ID","LMSD_HMDB_ID", "LMSD_SUB_CLASS",
                              "LMSD_CLASS_LEVEL4")
  
  # Create aggregate string of all lipid annotations for each LipidID
  sub_lmsd <- lmsd_ann_sub[, c("LipidID", "delta", "LMSD_NAME", "LMSD_SYSTEMATIC_NAME",
                               "LMSD_ABBREVIATION", "LMSD_SYNONYMS")]
  sub_lmsd$delta <- round(sub_lmsd$delta,6)
  agg_lmsd <- aggregate(. ~ LipidID, data = sub_lmsd, function(x) paste(unique(x), collapse = " ; "), na.action = na.pass)
  agg_lmsd[agg_lmsd == "NA"] <- NA
  rownames(agg_lmsd) <- agg_lmsd$LipidID
  agg_lmsd <- agg_lmsd[rownames(metab_SE),]
  
  # Select 1st lipid for multiple lipid matches (lowest delta should be at the top of list)
  top_LMSD_match <- lmsd_ann_sub %>%
    group_by(LipidID) %>%
    slice_head() %>%
    arrange(factor(LipidID, levels = rownames(Mass_data)))
  
  # Get the metabolite feature info data from SE object
  metab_info_temp <- data.frame(metab_SE@elementMetadata@listData, check.names = F, stringsAsFactors = T) %>% rownames_to_column(var = "LipidID")
  
  # Left join with to the main table using LipidID to get the correct ordering
  SE_metadata_added_lmsd <- metab_info_temp %>%
    left_join(top_LMSD_match, by = 'LipidID')
  rownames(SE_metadata_added_lmsd) <- SE_metadata_added_lmsd$LipidID
  
  lmsd_list <- list(
    "full_lmsd_ann" = distinct_lmsd_ann, ## full but with duplicates removed
    "agg_lmsd_df" = agg_lmsd,
    "metadata_lmsd_table" = SE_metadata_added_lmsd
  )
  
  return(lmsd_list)
}