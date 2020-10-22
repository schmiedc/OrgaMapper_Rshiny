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
                                   column_cell_table,
                                   column_orga_table) {
  
  # plots area, number of detections and mean value
  organelle_intensity_cell = ""
  
  if (plot_background_subtract) {
    
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
    
    if (plot_background_subtract) {
      
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

plot_detection_measurements <- function(data_table,
                                   column_cell_table,
                                   column_orga_table) {
  
}
