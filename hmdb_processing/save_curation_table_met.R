save_curation_table_met <- function(metab_SE, filename) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'readr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  df <- data.frame('Alignment_ID' = metab_SE@elementMetadata$`info.Alignment ID`,
                   'Metabolite_name' = metab_SE@elementMetadata$`info.Metabolite name`,
                   'Ionisation' = metab_SE@elementMetadata$ionisation,
                   'shortname' = metab_SE@elementMetadata$shortname,
                   'Fill' = metab_SE@elementMetadata$`info.Fill %`,
                   'S/N' = metab_SE@elementMetadata$`info.S/N average`,
                   'RT_avg' = metab_SE@elementMetadata$`info.Average Rt(min)`)
  write_csv(df, file = filename)
}