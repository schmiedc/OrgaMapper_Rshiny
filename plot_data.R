library(ggplot2)
# ==============================================================================
# plot detection and cell data
# ==============================================================================
boxplot_theme <- function() {
  theme_bw()
  
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # x axis ticks
    axis.text.x = element_text(color = "grey20", 
                               size = 10, 
                               angle = 45, 
                               hjust = .5, 
                               vjust = .5, 
                               face = "plain"),
    # y axis ticks
    axis.text.y = element_text(color = "grey20", 
                               size = 10, 
                               angle = 0, 
                               hjust = 1, 
                               vjust = 0, 
                               face = "plain"),
    # x axis labels
    axis.title.x = element_text(color = "grey20", 
                                size = 12, 
                                angle = 0, 
                                hjust = .5, 
                                vjust = 0, 
                                face = "plain"),
    # y axis labels
    axis.title.y = element_text(color = "grey20", 
                                size = 12, 
                                angle = 90, 
                                hjust = .5, 
                                vjust = .5, 
                                face = "plain"),
    # title
    title = element_text(color = "grey20", 
                         size = 14, 
                         angle = 0, 
                         hjust = 0, 
                         vjust = 1, 
                         face = "plain")
  )
}

lineplot_theme <- function() {
  theme_bw()
  
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # x axis ticks
    axis.text.x = element_text(color = "grey20", 
                               size = 10, 
                               angle = 0, 
                               hjust = .5, 
                               vjust = .5, 
                               face = "plain"),
    # y axis ticks
    axis.text.y = element_text(color = "grey20", 
                               size = 10, 
                               angle = 0, 
                               hjust = 1, 
                               vjust = 0, 
                               face = "plain"),
    # x axis labels
    axis.title.x = element_text(color = "grey20", 
                                size = 12, 
                                angle = 0, 
                                hjust = .5, 
                                vjust = 0, 
                                face = "plain"),
    # y axis labels
    axis.title.y = element_text(color = "grey20", 
                                size = 12, 
                                angle = 90, 
                                hjust = .5, 
                                vjust = .5, 
                                face = "plain"),
    # title
    title = element_text(color = "grey20", 
                         size = 14, 
                         angle = 0, 
                         hjust = 0, 
                         vjust = 1, 
                         face = "plain")
  )
}

plot_cell_measurements <- function(cell_data_table,
                                   plots,
                                   cell_column,
                                   orga_column,
                                   background_subtract) {
  
  # plots area, number of detections and mean value
  organelle_intensity_cell = ""
  
  if (background_subtract) {
    
    organelle_intensity_cell = "MeanOrgaBackSub"
    
  } else {
    
    organelle_intensity_cell = "MeanValueOrga"
    
  }
  
  measure1 <- list("Ferets", 
                   "CellArea", 
                   "NumDetections",
                   organelle_intensity_cell)
  
  measure1_title <- list("Average Feret's diameter", 
                         "Average cell area", 
                         "Number of detections per cell", 
                         "Average intensity in cell (detection channel)")
  
  measure1_label <- list("ferets diameter (µm)", 
                         "cell area (µm²)", 
                         "average count", 
                         "fluorescent intensity (A.U.)")
  
  measure1_file <- list("feretPlot", 
                        "cellArea", 
                        "numDetections", 
                        "intPerCellDetection")
  
  cell_measure_long <- cell_data_table %>% 
    pivot_longer(cols=Ferets:MeanOrgaBackSub,values_to = "measurement" )
  
  plot_list_cell <- list()
  
  for (index in seq_along(measure1)) {
    
    dataSubset <- subset(cell_measure_long, cell_measure_long$name==measure1[index])
    
    plot_cell <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
      geom_boxplot(outlier.size = 0, outlier.shape = 1) +
      stat_boxplot(geom = 'errorbar', width = 0.2) +
      geom_jitter(width = 0.1) +
      ggtitle(measure1_title[index]) + 
      xlab("Treatment") +
      ylab( measure1_label[index] ) + 
      boxplot_theme()
    
    plot_list_cell[[index]] <- plot_cell 
    
    ggsave(plot = plot_cell,
           file=paste0(plots, .Platform$file.sep, measure1_file[index], ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
  }
  
  # plot measure channel
  
  measure_intensity_cell = ""
  measure_intensity_peak = ""
  
  if (cell_column == 10 && orga_column == 10) {
    
    if (background_subtract) {
      
      measure_intensity_cell = "MeanMeasureBackSub"
      
    } else {
      
      measure_intensity_cell = "MeanValueMeasure"
      
    }

    cell_measure_filter_new <- cell_data_table[c("Name", measure_intensity_cell)]
    colnames(cell_measure_filter_new)[2] <- "measure"
    
    plot_cell_measure <- ggplot(cell_measure_filter_new, aes(x=Name, y=measure)) +
      ggtitle("Average intensity in cell (measure channel)") + 
      xlab("Treatment") +
      ylab("fluorescent intensity (A.U.)") +
      geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
      stat_boxplot(geom = 'errorbar', width = 0.2) +
      geom_jitter(width = 0.1) +
      boxplot_theme()
    
    plot_list_cell[[length(plot_list_cell)  + 1]] <- plot_cell_measure
    
    ggsave(plot = plot_cell_measure,
           file=paste0(plots, .Platform$file.sep, "intPerCellMeasure", ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
    
    }
    
    return (plot_list_cell)

}

plot_detection_measurements <- function(full_data_table,
                                        summary_table,
                                        plots,
                                        column_cell_table,
                                        column_orga_table,
                                        norm_distance_nucleus,
                                        background_subtract) {
  
  # lysosome density plots
  plot_list_detection <- list()
  
  name_count <- as.data.frame(table(full_data_table$Name))
  detect_list <- list()
  
  # goes through each experiment and calculates lysosome density
  # then peak normalizes the lysosome density
  # collects these normalized density plots in detect_list
  for (name in name_count$Var1){
    
    data_per_name <- subset(full_data_table, Name == name)
    density_per_name <- density(data_per_name$DistanceNorm, 
                                bw = "nrd0", 
                                n = 512, 
                                from = 0, 
                                to = norm_distance_nucleus)
    
    data_frame <- data.frame(density_per_name$x)
    data_frame$y <- density_per_name$y
    colnames(data_frame)[1] <-  "x"
    
    # peak normalisation
    max = max(data_frame$y, na.rm = FALSE)
    data_frame$peak_norm <- sapply(data_frame$y, function(x){x /  max})
    detect_list[[name]] <- data_frame
    
  }
  
  # binds collection of normalized density plots and binds them into one dataframe
  norm_list <- do.call("rbind", detect_list)
  norm_list1 <- tibble::rownames_to_column(norm_list, "nameindex")
  norm_list1_indices <- str_split_fixed(norm_list1$nameindex, "\\.", 2)
  norm_list2 <- cbind(norm_list1_indices, norm_list1)
  colnames(norm_list2)[1] <- "name"
  colnames(norm_list2)[2] <- "index"
  
  # Plot Lysosome density vs normalized distance from Nucleus
  # density plots without peak normalized data
  plot_density_raw <- ggplot(norm_list2, aes(x = x, 
                                             y = y, 
                                             group = name, 
                                             color = name)) + 
    geom_line() +
    xlab("Normalized distance from Nucleus") +
    ylab("Lysosome density (peak norm)") +
    ggtitle("Raw density plots") +
    scale_x_continuous(expand = c(0, 0)) + # force start at 0
    scale_y_continuous(expand = c(0, 0)) + # force start at 0
    lineplot_theme()
  
  ggsave(plot = plot_density_raw,
         file=paste0(plots, .Platform$file.sep, "raw_densityPlot", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_density_raw
  
  # Plot Lysosome density vs normalized distance from Nucleus
  # density plots with peak normalized data
  plot_density <- ggplot(norm_list2, aes(x = x, 
                                         y = peak_norm, 
                                         group = name, 
                                         color = name)) + 
    geom_line() +
    xlab("Normalized distance from Nucleus") +
    ylab("Lysosome density (peak norm)") +
    ggtitle("Peak normalized density plots") +
    scale_x_continuous(expand = c(0, 0)) + # force start at 0
    scale_y_continuous(expand = c(0, 0)) + # force start at 0
    lineplot_theme()
  
  ggsave(plot = plot_density,
         file=paste0(plots, .Platform$file.sep, "densityPlot", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_density
  
  # ------------------------------------------------------------------------------
  # plot distance and detection intensity
  organelle_intensity_peak = ""
  
  if (background_subtract) {
    
    organelle_intensity_peak = "PeakDetectBackSub.mean"
    
  } else {
    
    organelle_intensity_peak = "PeakDetectionInt.mean"
    
  }
  
  measure2 <- list("DistanceCal.mean", 
                   "DistanceNorm.mean",
                   organelle_intensity_peak)
  
  measure2_title <- list("Average distance from nucleus", 
                         "Average normalized distance from nucleus",
                         "Average peak detection intensity  (detection channel)")
  
  measure2_label <- list("distance (µm)", 
                         "normalized distance",
                         "fluorescent intensity (A.U.)")
  
  measure2_file <- list("distance", 
                        "distanceNorm",
                        "intPerDetectionDetectionChannel")
  
  summary_long <- summary_table %>% 
    pivot_longer(cols=DistanceRaw.mean:DistanceNorm.mean, 
                 values_to = "measurement" )
  
  
  for (index in seq_along(measure2)) {
    
    dataSubset <- subset(summary_long, summary_long$name==measure2[index])
    
    plot_detection <- distancePlot <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
      geom_boxplot(outlier.size = 0, outlier.shape = 1) +
      stat_boxplot(geom = 'errorbar', width = 0.2) +
      geom_jitter(width = 0.1) +
      ggtitle(measure2_title[index]) + 
      xlab("Treatment") +
      ylab(measure2_label[index]) +
      boxplot_theme()
    
    plot_list_detection[[length(plot_list_detection) + 1]] <- plot_detection
    
    ggsave(plot = plot_detection, 
           file=paste0(plots, .Platform$file.sep, measure2_file[index], ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
  }
  
  # ----------------------------------------------------------------------------
  # plot measure channel
  measure_intensity_peak = ""
  
  if (column_cell_table == 10 && column_orga_table == 10) {
    
    if (background_subtract) {
      
      measure_intensity_peak = "PeakMeasureBackSub.mean"
      
    } else {
      
      measure_intensity_peak = "PeakMeasureInt.mean"
      
    }
    
    summary_table_new <- summary_table[c("Name", measure_intensity_peak)]
    colnames(summary_table_new)[2] <- "measure"
    
    plot_peak_measure <- ggplot(summary_table_new, aes(x=Name, y=measure)) +
      ggtitle("Average peak detection intensity (measure channel)") + 
      xlab("Treatment") +
      ylab("fluorescent intensity (A.U.)") +
      geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
      stat_boxplot(geom = 'errorbar', width = 0.2) +
      geom_jitter(width = 0.1) +
      boxplot_theme()
    
    plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_peak_measure
    
    ggsave(plot = plot_peak_measure,
           file=paste0(plots, 
                       .Platform$file.sep, 
                       "intPerDetectionMeasureChannel", 
                       ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
  }
  
  return (plot_list_detection)
  
}
