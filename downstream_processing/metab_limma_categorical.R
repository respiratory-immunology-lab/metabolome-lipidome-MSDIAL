############################################################################################
# Copyright (c) 2021 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

# Define function to handle categorical limma testing
metab_limma_categorical <- function(metab_SE, metadata_var, metadata_condition = NULL, model_matrix = NULL, 
                                    contrast_matrix = NULL, adjust_method = 'BH', rownames = NULL, adj_pval_threshold = 0.05, 
                                    logFC_threshold = 1, legend_metadata_string = NULL,
                                    volc_plot_title = NULL, volc_plot_subtitle = NULL,
                                    volc_plot_xlab = NULL, volc_plot_ylab = NULL) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'limma', 'ggplot2', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Filter assay data by the conditional statement if present
  if(!is.null(metadata_condition)) {
    metab_SE <- metab_SE[, metadata_condition]
  }
  
  # Prepare limma data.frame
  metab_limma_data <- data.frame(assay(metab_SE))
  
  # Use shortname as rownames unless rownames provided
  if (is.null(rownames)) {
    rownames(metab_limma_data) <- make.unique(rowData(metab_SE)$shortname)
  } else {
    rownames(metab_limma_data) <- rownames
  }
  
  # Ensure that NA values are actually NA, and not character type 'NA'
  if (!is.null(metadata_condition)) {
    test_var <- metab_SE@metadata$metadata[[metadata_var]][metadata_condition]
  } else {
    test_var <- as.vector(metab_SE@metadata$metadata[[metadata_var]])
  }
  
  ensure_NA <- function(vector) {
    v <- replace(vector, vector == 'NA', NA)
    v <- replace(v, v == 'TRUE', 'Yes')
    v <- replace(v, v == 'FALSE', 'No')
    return(v)
  }
  
  test_var <- factor(ensure_NA(test_var))
  
  # Filter limma data to remove NA values
  metab_limma_data <- metab_limma_data[,!is.na(test_var)]
  
  # Create design matrix if none provided
  if (is.null(model_matrix)) {
    metab_limma_design <- model.matrix(~ 0 + test_var)
    colnames(metab_limma_design) <- levels(test_var)
  } else {
    metab_limma_design <- model_matrix
  }
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(metab_limma_data, metab_limma_design)
  
  # Define function to handle creation of the contrasts matrix
  make_contrasts_vector <- function(design_matrix) {
    # Define levels and make new lists
    levels <- colnames(design_matrix)
    num_levels <- length(levels)
    ci <- list()
    cj <- list()
    
    # For loops for generate contrast statements
    for (i in 2:num_levels) {
      ci[length(ci) + 1] <- paste0(levels[1], '-', levels[i])
      for (j in (i + 1):num_levels) {
        cj[length(cj) + 1] <- paste0(levels[i], '-', levels[j])
      }
    }
    
    # Unlist elements and remove NA values, then remove the last element (contrasts itself)
    cx <- c(unlist(ci), unlist(cj)); cx <- cx[!str_detect(cx, 'NA')]
    cx <- cx[1:length(cx) - 1]
    
    # Generate contrasts matrix and return it
    cont_matrix <- makeContrasts(contrasts = cx, levels = levels)
    cont_matrix
  }
  
  # Create the contrasts matrix if none provided
  if (is.null(contrast_matrix)) {
    cont_matrix <- make_contrasts_vector(metab_limma_design)
  } else {
    cont_matrix <- contrast_matrix
  }
  
  # Bayes statistics of differential expression
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = FALSE)
  
  # Generate a list of the differentially abundant metabolites
  get_all_topTables <- function(fit2, top = TRUE) {
    # Prepare an empty list object
    all_topTables <- list()
    
    # Run through each of the coefficients and add topTable to the list
    for (coef in 1:dim(fit2$contrasts)[2]) {
      coef_name <- colnames(fit2$contrasts)[coef]
      metab_topTable <- topTable(fit2, number = dim(fit2)[1], adjust.method = adjust_method, coef = coef)
      if (top == TRUE) {
        metab_topTable <- metab_topTable %>%
          filter(adj.P.Val < adj_pval_threshold) %>%
          filter(abs(logFC) >= logFC_threshold)
      }
      all_topTables[[coef_name]] <- metab_topTable
    }
    
    # Return list object
    all_topTables
  }
  
  # Get both significant results and all results lists
  metab_limma_signif <- get_all_topTables(fit2, top = TRUE)
  metab_limma_all <- get_all_topTables(fit2, top = FALSE)
  
  # Generate a venn diagram of the results
  results <- decideTests(fit2)
  venn_diagram <- vennDiagram(results)
  
  # Add a direction column to the 'metab_limma_all' topTables
  limma_add_categorical_direction <- function(metab_limma_all) {
    for (i in 1:length(metab_limma_all)) {
      metab_limma_all[[i]] <- data.frame(metab_limma_all[[i]]) %>%
        mutate(direction = case_when(
          adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ 'Decreased',
          adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ 'Increased',
          adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'NS'
        ))
    }
    metab_limma_all
  }
  
  metab_limma_all <- limma_add_categorical_direction(metab_limma_all)
  
  # Assign test_string if provided
  if (!is.null(legend_metadata_string)) {
    test_string <- legend_metadata_string
  } else {
    test_string <- 'Test Variable'
  }
  
  # Create a vector of strings for color legend
  plot_color_names <- function(test_string, metab_limma_all) {
    vector = ''
    for (i in 1:length(metab_limma_all)) {
      contrast_name <- names(metab_limma_all[i])
      contrast_name <- gsub('-', ' vs. ', contrast_name)
      contrast_name <- paste0(test_string, ':\n', contrast_name)
      vector <- c(vector, contrast_name)
    }
    vector <- vector[2:length(vector)]
    vector
  }
  
  plot_color_labels <- plot_color_names(test_string, metab_limma_all)
  
  # Define function to plot multiple volcano plots
  limma_volcano_plots <- function(topTable, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels) {
    # Create a blank list to hold plots
    plots <- list()
    
    # Make volcano plot for each topTable
    for (i in 1:length(metab_limma_all)) {
      plot <- ggplot(metab_limma_all[[i]], aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = direction)) +
        geom_hline(yintercept = -log10(adj_pval_threshold), linetype = 2) +
        geom_vline(xintercept = -logFC_threshold, linetype = 2) +
        geom_vline(xintercept = logFC_threshold, linetype = 2) +
        scale_color_manual(values = c('Decreased' = 'blue', 'Increased' = 'red', 'NS' = 'grey70'), name = plot_color_labels[i]) +
        geom_text_repel(data = metab_limma_all[[i]][metab_limma_all[[i]]$adj.P.Val < adj_pval_threshold & 
                                                      abs(metab_limma_all[[i]]$logFC) >= logFC_threshold, ],
                        label = rownames(metab_limma_all[[i]][metab_limma_all[[i]]$adj.P.Val < adj_pval_threshold & 
                                                                abs(metab_limma_all[[i]]$logFC) >= logFC_threshold, ]),
                        size = 2) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = plot_xlab,
             y = plot_ylab)
      
      # Add plot to list
      var_name <- names(metab_limma_all[i])
      plots[[var_name]] <- plot
    }
    
    # Return plots
    plots
  }
  
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
  
  # Plot all volcano plots
  volcano_plots <- limma_volcano_plots(metab_limma_all, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels)
  
  # Make a list of all components above to return from function
  return_list <- list(input_data = metab_limma_data,
                      test_variable = test_var,
                      model_matrix = metab_limma_design,
                      contrast_matrix = cont_matrix,
                      limma_significant = metab_limma_signif,
                      limma_all = metab_limma_all,
                      volcano_plots = volcano_plots,
                      venn_diagram = venn_diagram,
                      limma_type = 'categorical')
  
  # Return the list
  return_list
}