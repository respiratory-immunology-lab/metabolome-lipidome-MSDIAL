pmp_metab_isNamed <- function(metab_SE) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  row_data <- rowData(metab_SE)
  
  table(row_data$`info.Metabolite name` != 'Unknown')
}
