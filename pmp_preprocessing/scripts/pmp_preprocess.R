pmp_preprocess <- function(pos_df, neg_df, metadata = NULL, samples_key = 'Sample', intens_cols = NULL, info_cols = NULL,
                           blankFC = 5, max_perc_mv = 0.8, missingPeaksFraction = 0.8, max_rsd = 25, 
                           mv_imp_rowmax = 0.7, mv_imp_colmax = 0.7, mv_imp_method = 'knn'){
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'ggplot2', 'pmp', 'SummarizedExperiment', 'S4Vectors',
            'ggsci', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  metab_pos <- pos_df
  metab_neg <- neg_df
  
  # Get info columns and intensity columns
  if (is.null(info_cols)) {
    info_cols <- c('Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
                   'Adduct type', 'Post curation result', 'Fill %', 'MS/MS assigned',
                   'Reference RT', 'Reference m/z', 'Formula', 'Ontology',
                   'INCHIKEY', 'SMILES', 'Annotation tag (VS1.0)', 'RT matched',
                   'm/z matched', 'MS/MS matched', 'Comment', 'Manually modified for quantification',
                   'Manually modified for annotation', 'Isotope tracking parent ID',
                   'Isotope tracking weight number', 'Total score', 'RT similarity',
                   'Dot product', 'Reverse dot product', 'Fragment presence %', 'S/N average',
                   'Spectrum reference file name', 'MS1 isotopic spectrum', 'MS/MS spectrum')
    metab_pos_info <- metab_pos[, colnames(metab_pos) %in% info_cols]
    metab_neg_info <- metab_neg[, colnames(metab_neg) %in% info_cols]
  } else {
    metab_pos_info <- metab_pos[, info_cols]
    metab_neg_info <- metab_neg[, info_cols]
  }
  
  # Get QC, blanks, and sample columns
  if (is.null(intens_cols)) {
    intens_pattern <- paste('Blank_', 'QC_', paste0(samples_key, '_'), sep = '|')
    intens_cols <- str_detect(colnames(metab_pos), intens_pattern)
  }
  
  # Create dataframes for the SummarizedExperiment object
  metab_pos_counts <- as.matrix(metab_pos[, intens_cols])
  metab_neg_counts <- as.matrix(metab_neg[, intens_cols])
  
  # Remove MS/MS samples (not acquired in the same way)
  metab_pos_counts <- metab_pos_counts[, !(colnames(metab_pos_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  metab_neg_counts <- metab_neg_counts[, !(colnames(metab_neg_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  
  # Rename the data to indicate ionisation mode
  metab_pos_rownames <- paste0(metab_pos_info$`Alignment ID`, '_pos')
  metab_neg_rownames <- paste0(metab_neg_info$`Alignment ID`, '_neg')
  
  rownames(metab_pos_counts) <- metab_pos_rownames
  rownames(metab_neg_counts) <- metab_neg_rownames
  rownames(metab_pos_info) <- metab_pos_rownames
  rownames(metab_neg_info) <- metab_neg_rownames
  
  # Merge the positive and negative ionisation modes
  metab_counts <- rbind(metab_pos_counts, metab_neg_counts)
  metab_info <- rbind(metab_pos_info, metab_neg_info)
  
  # Create class and group vectors
  metab_class <- substr(colnames(metab_counts), start = 1, stop = 2)
  metab_group <- substr(colnames(metab_counts), start = 1, stop = 2)
  
  # Alternate steps depending on whether metadata was provided
  if (is.null(metadata)) {
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     rowData = list(info = metab_info),
                                     colData = DataFrame(class = metab_class))
  } else {
    # Check that the metadata matches the samples
    sample_cols <- stringr::str_detect(colnames(metab_counts), paste0(samples_key))
    metadata <- metadata[colnames(metab_counts)[sample_cols], ]
    identical(rownames(metadata), colnames(metab_counts)[sample_cols])
    
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     metadata = list(metadata = metadata),
                                     rowData = list(info = metab_info),
                                     colData = DataFrame(class = metab_class))
  }
  
  print('SummarizedExperiment object created...')
  
  ###
  ### FILTERING AND NORMALISATION
  ###
  
  # Original number of features
  features0 <- dim(metab_SE)
  
  # Replace missing values with NA to be compatible with downstream filtering
  assay(metab_SE) <- replace(assay(metab_SE), assay(metab_SE) == 0, NA)
  
  print('Replaced missing values with NA...')
  
  # Filter peaks and samples based on blanks
  metab_filt <- filter_peaks_by_blank(df = metab_SE,
                                      fold_change = blankFC,
                                      classes = metab_SE$class,
                                      remove_samples = TRUE,
                                      remove_peaks = TRUE,
                                      blank_label = 'Bl')
  
  print('Filtered peaks and samples based on blanks...')
  
  # Number of features
  features1 <- dim(metab_filt)
  
  # Filter samples based on missing values
  metab_filt <- filter_samples_by_mv(df = metab_filt,
                                     max_perc_mv = max_perc_mv)
  
  # Number of features
  features2 <- dim(metab_filt)
  
  # Filter peaks based on missing values
  metab_filt <- filter_peaks_by_fraction(df = metab_filt,
                                         min_frac = missingPeaksFraction,
                                         classes = metab_filt$class,
                                         method = 'across')
  
  # Number of features
  features3 <- dim(metab_filt)
  
  print('Filtered peaks and samples based on missing values...')
  
  # Filter peaks based on the % variation in the QC
  metab_filt <- filter_peaks_by_rsd(df = metab_filt,
                                    max_rsd = max_rsd,
                                    classes = metab_filt$class,
                                    qc_label = 'QC')
  
  # Number of features
  features4 <- dim(metab_filt)
  
  print('Filtered peaks based on the percentage variance in QC samples...')
  
  # Data normalisation
  metab_norm <- pqn_normalisation(df = metab_filt,
                                  classes = metab_filt$class,
                                  qc_label = 'QC')
  
  print('Normalised data using PQN...')
  print('Imputing values now...')
  
  # Missing values imputation
  metab_imp <- mv_imputation(df = metab_norm,
                             rowmax = mv_imp_rowmax,
                             colmax = mv_imp_colmax,
                             method = mv_imp_method)
  
  print('Finished imputing values...')
  
  # Data scaling
  metab_glog <- glog_transformation(df = metab_imp,
                                    classes = metab_imp$class,
                                    qc_label = 'QC')
  
  print('Scaled data via glog algorithm...')
  
  opt_lambda <- processing_history(metab_glog)$glog_transformation$lambda_opt
  
  glog_plot <- glog_plot_optimised_lambda(df = metab_imp,
                                          optimised_lambda = opt_lambda,
                                          classes = metab_imp$class,
                                          qc_label = 'QC')
  
  # Perform PCA
  PCA <- prcomp(t(assay(metab_glog)), center = TRUE)
  varexp <- c(summary(PCA)$importance[2,1]*100, summary(PCA)$importance[2,2]*100)
  
  # Create dataset
  data_PCA <- cbind(data.frame(Samples = rownames(PCA$x),
                               PC1 = PCA$x[,1],
                               PC2 = PCA$x[,2]),
                    class = metab_glog$class)
  
  # Plot results
  PCA_plot <- ggplot(data_PCA, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = class, color = class)) +
    stat_ellipse(aes(fill = class), geom = 'polygon', type = 't', level = 0.9, alpha = 0.2) +
    labs(title = 'Metabolomics',
         x = paste0('PC1 ', round(varexp[1], 2), '%'),
         y = paste0('PC2 ', round(varexp[2], 2), '%')) +
    scale_fill_jama(name = 'Class') +
    scale_color_jama(name = 'Class')
  
  # Make filtering dimensions dataframe
  filtering_dims <- data.frame(rbind(features0, features1, features2, features3, features4))
  colnames(filtering_dims) <- c('Features', 'Samples')
  rownames(filtering_dims) <- c('Original', 'Blank_filtered', 'MV_sample_filtered', 
                                'MV_peak_filtered', 'QC_var_filtered')
  
  # Prepare elements in a list to return from function
  results <- list(imputed_results = metab_imp,
                  glog_results = metab_glog,
                  glog_plot = glog_plot,
                  PCA_plot = PCA_plot,
                  filtering_dimensions = filtering_dims)
  
  print('Done!')
  
  results
}
