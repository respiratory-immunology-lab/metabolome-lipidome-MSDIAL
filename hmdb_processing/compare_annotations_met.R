compare_annotations_met <- function(metab_SE) {
  # Load packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_gnps_hmdb <- data.frame('Alignment.ID' = rownames(metab_SE),
                                 'Retention.Time' = rowData(metab_SE)$`info.Average Rt(min)`,
                                 'Fill_percent' = rowData(metab_SE)$`info.Fill %`,
                                 'S/N_ratio' = rowData(metab_SE)$`info.S/N average`,
                                 'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`,
                                 #'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps,
                                 'HMDB_annotation' = rowData(metab_SE)$HMDB,
                                 'HMDB_accession' = rowData(metab_SE)$HMDB_accession,
                                 'KEGG_annotation' = rowData(metab_SE)$KEGG)
  msdial_gnps_hmdb <- msdial_gnps_hmdb %>%
    column_to_rownames(var = 'Alignment.ID') %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA),
           KEGG_annotation= replace(KEGG_annotation, KEGG_annotation == '', NA)) %>%
    filter(!is.na(MSDIAL_annotation) | #!is.na(GNPS_annotation) | 
             !is.na(HMDB_annotation) | !is.na(KEGG_annotation))
  
  # Return the data.frame
  msdial_gnps_hmdb
}