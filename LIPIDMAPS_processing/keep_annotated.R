keep_annotated <- function(metab_SE) {
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$HMDB) | rowData(metab_SE)$`info.Metabolite name` != 'Unknown' |
                          !is.na(rowData(metab_SE)$LMSD_NAME) | !is.na(rowData(metab_SE)$LMSD_SYSTEMATIC_NAME) |
                          !is.na(rowData(metab_SE)$LMSD_ABBREVIATION)), ]
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  # Add names to "shortname" by preference: LMSD > HMDB > MS-DIAL > GNPS
  for (n in 1:nrow(metab_SE)) {
    if (!is.na(rowData(metab_SE)$LMSD_NAME[n]) | !is.na(rowData(metab_SE)$LMSD_SYSTEMATIC_NAME[n]) | !is.na(rowData(metab_SE)$LMSD_ABBREVIATION[n])) {
      rowData(metab_SE)$shortname[n] <- ifelse(is.na(rowData(metab_SE)$LMSD_NAME[n]),
                                               ifelse(is.na(rowData(metab_SE)$LMSD_ABBREVIATION[n]),
                                                      rowData(metab_SE)$LMSD_SYSTEMATIC_NAME[n],
                                                        rowData(metab_SE)$LMSD_ABBREVIATION[n]),
                                               rowData(metab_SE)$LMSD_NAME[n])
    } else if (!is.na(rowData(metab_SE)$HMDB[n]) )  {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$HMDB[n]
    } else if  (rowData(metab_SE)$`info.Metabolite name`[n] != "Unknown" & !(str_detect(rowData(metab_SE)$`info.Metabolite name`[n], "RIKEN"))) {
      rowData(metab_SE)$shortname[n] <- as.character(rowData(metab_SE)$`info.Metabolite name`[n])
    } else { ## This is for the RIKEN match - which will be labeled NA and removed
      rowData(metab_SE)$shortname[n] <- NA
    }
  }

  #Remove any NA's
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$shortname)), ]

  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '',
                                      gsub('; LC-ESI-QQ; MS2; CE', '',
                                           rowData(metab_SE)$shortname))
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub(' \\| .*', '', rowData(metab_SE)$shortname)
  rowData(metab_SE)$shortname <- gsub('([^;]*);.*', '\\1', rowData(metab_SE)$shortname)
  # Add the ionisation mode, and then make unique
  rowData(metab_SE)$shortname <- make.unique(paste0(rowData(metab_SE)$shortname,
                                                    '|', rowData(metab_SE)$ionisation))
  # Return the SE object
  return(metab_SE)
}
