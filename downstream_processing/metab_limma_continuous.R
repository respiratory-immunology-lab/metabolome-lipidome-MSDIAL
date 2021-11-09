############################################################################################
# Copyright (c) 2021 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

# Define custom function to handle continuous limma testing
metab_limma_continuous <- function(metab_SE, metadata_var, model_matrix = NULL, rownames = NULL, 
                                   adj_pval_threshold = 0.05, logFC_threshold = 1,
                                   volc_plot_title = NULL, volc_plot_subtitle = NULL,
                                   volc_plot_xlab = NULL, volc_plot_ylab = NULL) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'limma', 'ggplot2', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Prepare data.frame for limma (and remove any columns containing 'QC')
  metab_limma_data <- data.frame(assay(metab_SE)[, !str_detect(colnames(metab_SE), 'QC')])
  
  # Set rownames to shortname unless provided another vector
  if (is.null(rownames)) {
    rownames(metab_limma_data) <- make.unique(rowData(metab_SE)$shortname)
  } else {
    rownames(metab_limma_data) <- rownames
  }
  
  # Retrieve the test variable from the metadata if given a string
  if (class(metadata_var) == 'character') {
    test_var <- metab_SE@metadata$metadata[[metadata_var]]
  } else {
    test_var <- metadata_var # otherwise use metadata_var directly
  }
  
  # Create design matrix unless provided one
  if (is.null(model_matrix)) {
    metab_limma_design <- model.matrix(~ test_var)
  } else {
    metab_limma_design
  }
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(metab_limma_data, metab_limma_design)
  
  # Bayes statistics of differential expression
  fit <- eBayes(fit, robust = TRUE, trend = FALSE)
  
  # Generate a list of the top differentially abundant metabolites
  metab_limma_top <- topTable(fit, number = dim(fit)[1], adjust = 'BH') %>%
    filter(adj.P.Val < adj_pval_threshold) %>%
    filter(abs(logFC) >= logFC_threshold)
  
  # Generate a list of all differentially abundant metabolites
  metab_limma_all <- topTable(fit, number = dim(fit)[1], adjust = 'BH') %>%
    mutate(direction = case_when(
      adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ 'Decreased',
      adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ 'Increased',
      adj.P.Val >= adj_pval_threshold | abs(logFC) < 1 ~ 'NS'
    ))
  
  # Prepare labelling
  if (is.null(volc_plot_title)) {
    plot_title <- 'Differential Intensity Metabolites'
  } else {
    plot_title <- volc_plot_title
  }
  
  if (is.null(volc_plot_subtitle)) {
    plot_subtitle <- NULL
  } else {
    plot_subtitle <- volc_plot_subtitle
  }
  
  if (is.null(volc_plot_xlab)) {
    plot_xlab <- 'log2FC'
  } else {
    plot_xlab <- volc_plot_xlab
  }
  
  if (is.null(volc_plot_ylab)) {
    plot_ylab <- 'Significance\n-log10(adjusted p-value)'
  } else {
    plot_ylab <- volc_plot_ylab
  }
  
  # Plot volcano
  metab_limma_volc <- ggplot(metab_limma_all, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = direction)) +
      geom_hline(yintercept = -log10(adj_pval_threshold), linetype = 2) +
      geom_vline(xintercept = -logFC_threshold, linetype = 2) +
      geom_vline(xintercept = logFC_threshold, linetype = 2) +
      scale_color_manual(values = c('Decreased' = 'blue', 'Increased' = 'red', 'NS' = 'grey70'), name = 'Direction') +
      geom_text_repel(data = metab_limma_all[metab_limma_all$adj.P.Val < adj_pval_threshold & abs(metab_limma_all$logFC) >= logFC_threshold, ],
                      label = rownames(metab_limma_all[metab_limma_all$adj.P.Val < adj_pval_threshold & abs(metab_limma_all$logFC) >= logFC_threshold, ]),
                      size = 2) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = plot_xlab,
           y = plot_ylab)
  
  # Prepare a list object to return various elements
  return_list <- list(input_data = metab_limma_data,
                      test_variable = test_var,
                      limma_significant = metab_limma_top,
                      all_values = metab_limma_all,
                      volcano_plot = metab_limma_volc,
                      limma_type = 'continuous')
  
  # Return the list
  return_list
}