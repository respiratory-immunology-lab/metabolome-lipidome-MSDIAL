save_curation_table <- function(metab_SE, filename) {
    # Load packages
    pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
    for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }

    # Prepare a data.frame with useful information for curation
    df <- data.frame('Alignment_ID' = rowData(metab_SE)$`info.Alignment ID`,
                    'Fill_percent' = rowData(metab_SE)$`info.Fill %`,
                    'S/N_average' = rowData(metab_SE)$`info.S/N average`,
                    'Metabolite_name' = rowData(metab_SE)$`info.Metabolite name`,
                    'Ionisation' = rowData(metab_SE)$ionisation,
                    'shortname' = rowData(metab_SE)$shortname,
                    'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps,
                    'HMDB_annotation' = rowData(metab_SE)$HMDB,
                    'HMDB_accession' = rowData(metab_SE)$HMDB_accession,
                    'KEGG_annotation' = rowData(metab_SE)$KEGG)

    # Write the data.frame to a .csv file at the requested location
    write_csv(df, file = filename)
}