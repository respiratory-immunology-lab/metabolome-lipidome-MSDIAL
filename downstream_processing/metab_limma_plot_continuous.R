############################################################################################
# Copyright (c) 2021 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for more information)     #
############################################################################################

# Plot continuous metadata variable vs metabolite
metab_limma_plot_indiv_continuous <- function(metab_limma_cont_object, feature_num, plot_subtitle = NULL, 
                                              plot_x_label = NULL, plot_y_label = NULL, 
                                              text_size = 8, geom_point_fill = 'black', geom_point_alpha = 0.7, geom_point_size = 2) {
  shortname <- rownames(metab_limma_cont_object$limma_significant)[feature_num]
  
  logFC <- metab_limma_cont_object$limma_significant$logFC[feature_num]
  direction <- ifelse(logFC > 0, 'Increased', 'Decreased')
  
  metab_df <- data.frame(t(metab_limma_cont_object$input_data[shortname,]),
                         metab_limma_cont_object$test_variable[!is.na(metab_limma_cont_object$test_variable)],
                         direction)
  
  colnames(metab_df) <- c('metab', 'var', 'direction')
  
  if (is.null(plot_y_label)) {
    plot_y_label <- 'Intensity'
  }
  
  plot <- ggplot(metab_df, aes(x = var, y = metab)) +
    geom_point(shape = 21, fill = 'black', alpha = 0.7, size = geom_point_size) +
    geom_smooth(aes(color = direction), method = 'gam', formula = y ~ s(x, bs = 'cs')) +
    labs(title = shortname,
         subtitle = plot_subtitle,
         x = plot_x_label,
         y = plot_y_label) +
    theme(legend.position = 'NONE',
          text = element_text(size = text_size)) +
    scale_color_manual(values = c('Increased' = 'red', 'Decreased' = 'blue'))
  
  return(plot)
}

# Plot continuous metadata variable vs all metabolites
metab_limma_plot_all_continuous <- function(metab_limma_cont_object, plot_subtitle = NULL, 
                                            plot_x_label = NULL, plot_y_label = NULL, 
                                            text_size = 8, geom_point_fill = NULL, geom_point_alpha = NULL,
                                            geom_point_size = 2) {
  # Create an empty list to hold the plots
  plot_list <- list()
  
  # Create plots for each significant metabolite
  for (metab in 1:dim(metab_limma_cont_object$limma_significant)[1]) {
    shortname <- rownames(metab_limma_cont_object$limma_significant)[metab]
    
    p <- metab_limma_plot_indiv_continuous(metab_limma_cont_object = metab_limma_cont_object,
                                           feature_num = metab,
                                           plot_subtitle = plot_subtitle,
                                           plot_x_label = plot_x_label,
                                           plot_y_label = plot_y_label,
                                           text_size = text_size,
                                           geom_point_fill = geom_point_fill,
                                           geom_point_alpha = geom_point_alpha,
                                           geom_point_size = geom_point_size)
    
    plot_list[[shortname]] <- p
  }
  
  return(plot_list)
}