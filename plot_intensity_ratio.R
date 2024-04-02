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

# compute intensity ratio
# needs to be revised
compute_intensity_ratio <- function(dataframe_value, perimeter, bin_width, perimeter_offset) {

  
  perimeter_periphery = ( perimeter + perimeter_offset ) / bin_width
  perimeter_bin = perimeter / bin_width

  value_perinuclear <- dataframe_value %>% 
    filter(intensityDistanceCalibrated <= perimeter_bin) %>% 
    group_by(identifier, series, cell) %>%  
    summarise(mean_perinuclear = mean(orgaIntensityBacksub), na.rm =TRUE)
  
  value_peripheral <- dataframe_value %>% 
    drop_na() %>% 
    filter(intensityDistanceCalibrated > perimeter_periphery) %>% 
    group_by(identifier, series, cell) %>% 
    summarise(mean_peripheral = mean(orgaIntensityBacksub))
  
  value_merge <- inner_join(value_perinuclear, value_peripheral, by = c("identifier", "series", "cell"))
  value_merge$intensity_ratio <- value_merge$mean_perinuclear / value_merge$mean_peripheral
  
  return (value_merge)
  
}

plot_intensity_ratio <- function(dataFrame, name, directory) {
  
  plot_intensity <- ggplot(dataFrame, aes(x=identifier, y=intensity_ratio)) +
    ggtitle(sprintf("Intensity Ratio \n%s channel", name)) + 
    xlab("Treatment") +
    ylab("I perinculear / I peripheral") +
    geom_boxplot(outlier.size = 0, outlier.shape = NA, na.rm=TRUE) + 
    stat_boxplot(geom = 'errorbar', width = 0.2, na.rm=TRUE) +
    # scale_y_continuous(limits = quantile(dataFrame$intensity_ratio, c(0.1, 0.9))) +
    geom_jitter(width = 0.1, na.rm=TRUE) + 
    boxplot_theme()
  
  ggsave(plot = plot_intensity,
         file=paste0(directory, .Platform$file.sep, sprintf("%s_intensityRatio", name), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")

}