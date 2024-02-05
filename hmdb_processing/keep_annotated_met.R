keep_annotated_met <- function(metab_SE) {
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$HMDB) | #!is.na(rowData(metab_SE)$compound_name_gnps) |
                          rowData(metab_SE)$`info.Metabolite name` != 'Unknown'), ]
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  
  # Add names to "shortname" by preference: HMDB > GNPS > MS-DIAL
  for (n in 1:nrow(metab_SE)) {
    if (!is.na(rowData(metab_SE)$HMDB)[n]) {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$HMDB[n]
    } #else if (!is.na(rowData(metab_SE)$compound_name_gnps)[n]) {
    #rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$compound_name_gnps[n]
    else {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$`info.Metabolite name`[n]
    }
  }
  
  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '', 
                                      gsub('; LC-ESI-QQ; MS2; CE', '', 
                                           rowData(metab_SE)$shortname))
  # Remove leftover trailing whitespace
  rowData(metab_SE)$shortname <- trimws(rowData(metab_SE)$shortname)
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub('([^;]*);.*', '\\1', rowData(metab_SE)$shortname)
  # Add the ionisation mode, and then make unique
  rowData(metab_SE)$shortname <- make.unique(paste0(rowData(metab_SE)$shortname, 
                                                    '_', rowData(metab_SE)$ionisation))
  # Return the SE object
  metab_SE
}