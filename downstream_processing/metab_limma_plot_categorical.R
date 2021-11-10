############################################################################################
# Copyright (c) 2021 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

# Define function to plot individual significant features from the limma categorical wrapper function
metab_limma_plot_indiv_categorical <- function(metab_limma_cat_object, metab_limma_cat_comparison, feature_num, 
                                               plot_subtitle = NULL, plot_x_lab = NULL, plot_y_lab = NULL, 
                                               plot_stat_comparisons, plot_fill_name = NULL,
                                               text_size = 8, pval_text_size = 2) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'ggplot2', 'ggpubr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }                                               

  shortname <- rownames(metab_limma_cat_object$limma_significant[[metab_limma_cat_comparison]])[feature_num]
  
  metab_df <- data.frame(t(metab_limma_cat_object$input_data[shortname,]),
                         metab_limma_cat_object$test_variable[!is.na(metab_limma_cat_object$test_variable)])
  
  colnames(metab_df) <- c('metab', 'var')
  
  plot <- ggplot(metab_df, aes(x = var, y = metab)) +
    geom_boxplot(aes(y = metab, fill = var), width = 0.6, alpha = 0.2) +
    geom_dotplot(aes(fill = var), binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
    labs(title = shortname,
         subtitle = plot_subtitle,
         x = plot_x_lab,
         y = plot_y_lab) +
    scale_fill_jama(name = plot_fill_name) +
    theme(legend.position = 'NONE',
          text = element_text(size = text_size)) +
    stat_compare_means(comparisons = plot_stat_comparisons, size = pval_text_size)
  
  return(plot)
}

# Define function to make plots for all significant features from the limma categorical wrapper function
metab_limma_plot_all_categorical <- function(metab_limma_cat_object, metab_limma_cat_comparison, plot_subtitle = NULL, 
                                            plot_x_lab = NULL, plot_y_lab = NULL, plot_stat_comparisons, plot_fill_name = NULL,
                                            text_size = 8, pval_text_size = 2) {
  # Create an empty list to hold the plots
  plot_list <- list()
  
  for (i in 1:dim(metab_limma_cat_object$limma_significant[[metab_limma_cat_comparison]])[1]) {
    # Generate the individual plot
    p <- metab_limma_plot_indiv_categorical(metab_limma_cat_object = metab_limma_cat_object,
                                            metab_limma_cat_comparison = metab_limma_cat_comparison,
                                            feature_num = i,
                                            plot_subtitle = plot_subtitle,
                                            plot_x_lab = plot_x_lab,
                                            plot_y_lab = plot_y_lab,
                                            plot_stat_comparisons = plot_stat_comparisons,
                                            plot_fill_name = plot_fill_name,
                                            text_size = text_size,
                                            pval_text_size = pval_text_size)
    # Add the plot to the plot list
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}