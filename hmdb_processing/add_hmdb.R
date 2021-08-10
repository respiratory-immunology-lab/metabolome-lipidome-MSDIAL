# Function to add HMDB information
add_hmdb <- function(metab_SE, hmdb, mass_tol = 0.002) {
  # Set ion mass
  ion <- 1.007276
  # Transform everything into a vector for faster looping
  hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
  hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  Mass_db <- as.vector(as.numeric(hmdb$monisotopic_molecular_weight))
  KEGG_db <- as.vector(hmdb$kegg_id)
  Name_db <- as.vector(hmdb$name)
  Mass_data <- as.vector(rowData(metab_SE)$`info.Average Mz`)
  Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  HMDB_data <- rep(NA, length(Mass_data))
  KEGG_data <- rep(NA, length(Mass_data))
  # Get masses corrected for ion precursors
  for (m in 1:length(Mass_data)){
    if (Adduct_data[m] %in% c('[M+2H]2+','[M+H]+','[M+Na]+','[M+NH4]+','[M+NH4]2+')){
      Mass_data[m] <- Mass_data[m]-ion}
    else {
      Mass_data[m] <- Mass_data[m]+ion}}
  # Run loop
  for (n in 1:length(Mass_db)){
    for (m in 1:length(Mass_data)){
      if (is.na(HMDB_data[m])==TRUE & between(Mass_db[n], Mass_data[m]-mass_tol, Mass_data[m]+mass_tol)==TRUE){
        HMDB_data[m] <- Name_db[n]
        KEGG_data[m] <- KEGG_db[n]}
      else if (is.na(HMDB_data[m])==FALSE & between(Mass_db[n], Mass_data[m]-mass_tol, Mass_data[m]+mass_tol)==TRUE){
        HMDB_data[m] <- paste(HMDB_data[m],Name_db[n], sep=';')
        KEGG_data[m] <- paste(KEGG_data[m],KEGG_db[n], sep=';')}}}
  # Add new information to SE experiment object
  rowData(metab_SE)$HMDB <- HMDB_data
  rowData(metab_SE)$KEGG <- KEGG_data
  return(metab_SE)}
