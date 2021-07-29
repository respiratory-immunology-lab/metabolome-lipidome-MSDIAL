gnps_SE_names <- function(gnps_df, metab_SE) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Get the metabolite feature info data from SE object
  metab_info_temp <- as.data.frame(metab_SE@elementMetadata@listData)
  metab_info_temp$alignment_ionisation <- rownames(metab_info_temp)
  
  # Select just the GNPS name (and alignment_ionisation column for joining)
  gnps_to_match <- gnps_df %>% dplyr::select(alignment_ionisation, compound_name_gnps)
  
  # Left join the GNPS names to the main table in order to get the right order
  metab_info_temp <- metab_info_temp %>% 
    left_join(gnps_to_match, by = 'alignment_ionisation')
  
  # Obtain the vector with GNPS names
  compound_name_gnps <- metab_info_temp$compound_name_gnps
  
  # Return the vector
  compound_name_gnps
}
