compare_annotations_lip <- function(metab_SE, agg_lmsd_ann) {
  # Load packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # LMSD_annotation will prioritise use of LMSD_NAME,then LMSD_ABBREVIATION and lastly LMSD_SYSTEMATIC_NAME
  LMSD_ann <- data.frame("ann" = ifelse(is.na(agg_lmsd_ann$LMSD_NAME),
                                        ifelse(is.na(agg_lmsd_ann$LMSD_ABBREVIATION),
                                               ifelse(is.na(agg_lmsd_ann$LMSD_SYSTEMATIC_NAME), NA, agg_lmsd_ann$LMSD_SYSTEMATIC_NAME),
                                               agg_lmsd_ann$LMSD_ABBREVIATION),
                                        agg_lmsd_ann$LMSD_NAME))
  rownames(LMSD_ann) <- agg_lmsd_ann$LipidID
  
  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_lmsd_gnps_hmdb <- data.frame('LipidID' = rowData(metab_SE)$LipidID,
                                      'Mz' = rowData(metab_SE)$`info.Average Mz`,
                                      'RT' = rowData(metab_SE)$`info.Average Rt(min)`,
                                      'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`,
                                      #'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps,
                                      #'HMDB_annotation' = rowData(metab_SE)$HMDB,
                                      #'KEGG_annotation' = rowData(metab_SE)$KEGG
                                      'LMSD_annotation' = LMSD_ann$ann
  ) %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA) #, 
           #KEGG_annotation= replace(KEGG_annotation, KEGG_annotation == '', NA)
    ) %>%
    filter(!is.na(MSDIAL_annotation) | !is.na(LMSD_annotation)
           #!is.na(GNPS_annotation) | 
           #!is.na(HMDB_annotation) | 
           #!is.na(KEGG_annotation) | 
    )
  
  # Return the data.frame
  return(msdial_lmsd_gnps_hmdb)
}