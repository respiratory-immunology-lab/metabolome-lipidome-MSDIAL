############################################################################################
# Copyright (c) 2021 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

### HEATMAP FOR METAB_LIMMA_CONTINUOUS OUTPUT ###
metab_limma_plot_heatmap_continuous <- function(metab_limma_cont_object, 
                                                metadata_to_include = NULL, 
                                                column_annotation_labels = NULL,
                                                continuous_colour_ramp = NULL,
                                                heatmap_scale_name = 'Intensity',
                                                heatmap_column_title = 'Differential Intensity Metabolites',
                                                heatmap_row_title = 'Feature',
                                                heatmap_order_columns_by_test_variable = TRUE,
                                                heatmap_show_column_names = TRUE,
                                                heatmap_rowname_text_size = 8,
                                                heatmap_annotation_legend_param = NULL,
                                                save_to_pdf = TRUE,
                                                save_to_png = TRUE,
                                                output_filename = NULL) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'grid', 'ComplexHeatmap', 'circlize')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Retrieve the test variable
  test_variable <- metab_limma_cont_object$test_variable
  
  # Retrieve additional metadata values if not null
  if (!is.null(metadata_to_include)) {
    metadata_variables <- metab_limma_cont_object$input_metadata %>%
      dplyr::select(all_of(metadata_to_include))
  }
  
  # Create a data.frame for the column annotation
  if (!is.null(metadata_to_include)) {
    colAnn <- data.frame(cbind(metadata_variables, test_variable))
  } else {
    colAnn <- data.frame(test_variable)
  }
  if (!is.null(column_annotation_labels) & length(column_annotation_labels) == ncol(colAnn)) {
    colnames(colAnn) <- column_annotation_labels
  }
  
  # Define colours
  if (!is.null(continuous_colour_ramp) & class(continuous_colour_ramp) == 'function') {
    continuous_colours <- continuous_colour_ramp
  } else {
    continuous_colours <- circlize::colorRamp2(c(0, max(test_variable) / 2, max(test_variable)), c('goldenrod', 'cadetblue1', 'dodgerblue4'))
  }

  last_col <- colnames(colAnn)[length(colnames(colAnn))] # the test variable is in the last column
  colours <- list() # create an empty list
  colours[[last_col]] <- continuous_colours # use the name of the test variable as the list item name
  
  # Make heatmap annotation object
  if (!is.null(heatmap_annotation_legend_param)) {
    colAnn_hm <- HeatmapAnnotation(df = colAnn, col = colours, annotation_legend_param = heatmap_annotation_legend_param)
  } else {
    colAnn_hm <- HeatmapAnnotation(df = colAnn, col = colours)
  }
  
  # Retrieve names of significant features
  cont_sig_features <- rownames(metab_limma_cont_object$limma_significant)  
    
  # Subset input data for significant features only
  input_data_subset <- metab_limma_cont_object$input_data[cont_sig_features, ]
  
  # Set heatmap parameters
  if (heatmap_order_columns_by_test_variable == TRUE) {
    column_order = order(metab_limma_cont_object$test_variable)
  } else {
    column_order = NULL
  }
  
  # Make the heatmap
  limma_cont_heatmap <- ComplexHeatmap::Heatmap(matrix = t(scale(t(input_data_subset))),
                                                name = heatmap_scale_name,
                                                column_title = heatmap_column_title,
                                                row_title = heatmap_row_title,
                                                column_labels = colnames(input_data_subset),
                                                top_annotation = colAnn_hm,
                                                rect_gp = gpar(col = 'white', lwd = 0.9),
                                                column_order = column_order,
                                                show_column_names = heatmap_show_column_names,
                                                row_names_gp = gpar(fontsize = heatmap_rowname_text_size))
  
  # Save to files if requested
  if (isTRUE(save_to_pdf) & !is.null(output_filename)) {
    pdf(paste0(output_filename, '.pdf'), 
        width = 0.1*ncol(input_data_subset) + 4, 
        height = 0.1*nrow(input_data_subset) + 1)
    draw(limma_cont_heatmap)
    dev.off()
  }
  if (isTRUE(save_to_png) & !is.null(output_filename)) {
    png(paste0(output_filename, '.png'), 
        width = 0.1*ncol(input_data_subset) + 4, 
        height = 0.1*nrow(input_data_subset) + 1,
        units = 'in', res = 300)
    draw(limma_cont_heatmap)
    dev.off()
  }
  
  # Prepare a grob version of the heatmap
  grob_hm <- grid::grid.grabExpr(draw(limma_cont_heatmap))
  
  # Prepare return list
  return_list <- list(heatmap = limma_cont_heatmap,
                      heatmap_grob = grob_hm)
  
  return(return_list)
}

### HEATMAP FOR METAB_LIMMA_CATEGORICAL OUTPUT ###
metab_limma_plot_heatmap_categorical <- function(metab_limma_cat_object, 
                                                 metadata_to_include = NULL, 
                                                 column_annotation_labels = NULL,
                                                 categorical_colours = NULL, # must be a named vector of length equal to test factor levels
                                                 heatmap_scale_name = 'Intensity',
                                                 heatmap_column_title = 'Differential Intensity Metabolites',
                                                 heatmap_row_title = 'Feature',
                                                 heatmap_order_columns_by_test_variable = TRUE,
                                                 heatmap_show_column_names = TRUE,
                                                 heatmap_rowname_text_size = 8,
                                                 heatmap_annotation_legend_param = NULL,
                                                 save_to_pdf = TRUE,
                                                 save_to_png = TRUE,
                                                 output_filename = NULL) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'grid', 'ComplexHeatmap', 'circlize')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Retrieve the test variable
  test_variable <- metab_limma_cat_object$test_variable[!is.na(metab_limma_cat_object$test_variable)]
  
  # Retrieve additional metadata values if not null
  if (!is.null(metadata_to_include)) {
    metadata_variables <- metab_limma_cat_object$input_metadata %>%
      dplyr::select(all_of(metadata_to_include))
  }
  
  # Create a data.frame for the column annotation
  if (!is.null(metadata_to_include)) {
    colAnn <- data.frame(cbind(metadata_variables, test_variable))
  } else {
    colAnn <- data.frame(test_variable)
  }
  if (!is.null(column_annotation_labels) & length(column_annotation_labels) == ncol(colAnn)) {
    colnames(colAnn) <- column_annotation_labels
  }
  
  # Define colours
  if (!is.null(categorical_colours) & length(categorical_colours) == length(levels(test_variable))) {
    categorical_colours <- categorical_colours
  }

  last_col <- colnames(colAnn)[length(colnames(colAnn))] # the test variable is in the last column
  colours <- list() # create an empty list
  colours[[last_col]] <- categorical_colours # use the name of the test variable as the list item name
  
  # Make heatmap annotation object
  if (!is.null(heatmap_annotation_legend_param)) {
    colAnn_hm <- HeatmapAnnotation(df = colAnn, col = colours, annotation_legend_param = heatmap_annotation_legend_param)
  } else {
    colAnn_hm <- HeatmapAnnotation(df = colAnn, col = colours)
  }
  
  # Retrieve the categorical comparisons
  cat_comparisons <- names(metab_limma_cat_object$limma_significant)
  
  # Create empty list for heatmaps
  heatmap_list <- list()
  
  # Loop through the categorical comparisons and make heatmaps for each
  for (comp in cat_comparisons) {
    # Retrieve names of significant features
    cat_sig_features <- rownames(metab_limma_cat_object$limma_significant[[comp]])  
      
    # Subset input data for significant features only
    input_data_subset <- metab_limma_cat_object$input_data[cat_sig_features, ]
    
    # Continue only if there are actually rows present
    if (nrow(input_data_subset) > 0) {
      # Set heatmap parameters
      if (heatmap_order_columns_by_test_variable == TRUE) {
        column_order = order(metab_limma_cat_object$test_variable[!is.na(metab_limma_cat_object$test_variable)])
      } else {
        column_order = NULL
      }
      
      # Make the heatmap
      limma_cat_heatmap <- ComplexHeatmap::Heatmap(matrix = t(scale(t(input_data_subset))),
                                                    name = heatmap_scale_name,
                                                    column_title = heatmap_column_title,
                                                    row_title = heatmap_row_title,
                                                    column_labels = colnames(input_data_subset),
                                                    top_annotation = colAnn_hm,
                                                    rect_gp = gpar(col = 'white', lwd = 0.9),
                                                    column_order = column_order,
                                                    show_column_names = heatmap_show_column_names,
                                                    row_names_gp = gpar(fontsize = heatmap_rowname_text_size))
      
      # Add heatmap to list
      heatmap_list[[comp]] <- limma_cat_heatmap
      
      # Save to files if requested
      if (isTRUE(save_to_pdf) & !is.null(output_filename)) {
        pdf(paste0(output_filename, comp, '.pdf'), 
            width = 0.1*ncol(input_data_subset) + 4, 
            height = 0.1*nrow(input_data_subset) + 1)
        draw(limma_cat_heatmap)
        dev.off()
      }
      if (isTRUE(save_to_png) & !is.null(output_filename)) {
        png(paste0(output_filename, comp, '.png'), 
            width = 0.1*ncol(input_data_subset) + 4, 
            height = 0.1*nrow(input_data_subset) + 1,
            units = 'in', res = 300)
        draw(limma_cat_heatmap)
        dev.off()
      }
    }
  }
  
  # Prepare grob versions of the heatmaps
  grob_hms <- list()
  for (hm in names(heatmap_list)) {
    grob_hms[[hm]] <- grid::grid.grabExpr(draw(heatmap_list[[hm]]))
  }
  
  # Prepare return list
  return_list <- list(heatmaps = heatmap_list,
                      heatmap_grobs = grob_hms)
  
  return(return_list)
}

### OUTER FUNCTION ###
metab_limma_plot_heatmap <- function(metab_limma_object, 
                                     metadata_to_include = NULL, 
                                     column_annotation_labels = NULL,
                                     continuous_colour_ramp = NULL, # made using the circlize::colorRamp2 function
                                     categorical_colours = NULL, # must be a named vector of length equal to test factor levels
                                     heatmap_scale_name = 'Intensity',
                                     heatmap_column_title = 'Differential Intensity Metabolites',
                                     heatmap_row_title = 'Feature',
                                     heatmap_order_columns_by_test_variable = TRUE,
                                     heatmap_show_column_names = TRUE,
                                     heatmap_rowname_text_size = 8,
                                     heatmap_annotation_legend_param = NULL,
                                     save_to_pdf = TRUE,
                                     save_to_png = TRUE,
                                     output_filename = NULL) {
  # Determine which inner function to use based on limma_type
  if (metab_limma_object$limma_type == 'continuous') {
    output <- metab_limma_plot_heatmap_continuous(metab_limma_cont_object = metab_limma_object, 
                                                  metadata_to_include = metadata_to_include, 
                                                  column_annotation_labels = column_annotation_labels,
                                                  continuous_colour_ramp = continuous_colour_ramp,
                                                  heatmap_scale_name = heatmap_scale_name,
                                                  heatmap_column_title = heatmap_column_title,
                                                  heatmap_row_title = heatmap_row_title,
                                                  heatmap_order_columns_by_test_variable = heatmap_order_columns_by_test_variable,
                                                  heatmap_show_column_names = heatmap_show_column_names,
                                                  heatmap_rowname_text_size = heatmap_rowname_text_size,
                                                  heatmap_annotation_legend_param = heatmap_annotation_legend_param,
                                                  save_to_pdf = save_to_pdf,
                                                  save_to_png = save_to_png,
                                                  output_filename = output_filename)
  }
  if (metab_limma_object$limma_type == 'categorical') {
    output <- metab_limma_plot_heatmap_categorical(metab_limma_cat_object = metab_limma_object, 
                                                   metadata_to_include = metadata_to_include, 
                                                   column_annotation_labels = column_annotation_labels,
                                                   categorical_colours = categorical_colours, # must be a named vector of length equal to test factor levels
                                                   heatmap_scale_name = heatmap_scale_name,
                                                   heatmap_column_title = heatmap_column_title,
                                                   heatmap_row_title = heatmap_row_title,
                                                   heatmap_order_columns_by_test_variable = heatmap_order_columns_by_test_variable,
                                                   heatmap_show_column_names = heatmap_show_column_names,
                                                   heatmap_rowname_text_size = heatmap_rowname_text_size,
                                                   heatmap_annotation_legend_param = heatmap_annotation_legend_param,
                                                   save_to_pdf = save_to_pdf,
                                                   save_to_png = save_to_png,
                                                   output_filename = output_filename)
  }
  return(output)
}