library(ggplot2)
library(tidyverse)

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
# ------------------------------------------------------------------------------
# create binned data 
plot_intensity_map <- function(dataframe_calib, 
                               dataframe_norm, 
                               channel_string, 
                               bin_width_calib, 
                               upper_limit_calib, 
                               bin_width_normalized,
                               upper_limit_normalized,
                               plots) {
  
  if (channel_string == "orga") {
    
    bin_param_calib = "orga_mean"
    bin_param_norm = "orga_mean"
    
  } else if (channel_string == "measure") {
    
    bin_param_calib = "measure_mean"
    bin_param_norm = "measure_mean"
    
  }
  
  identifier_count <- as.data.frame(table(dataframe_calib$identifier))
  
  subset_list <- list()
  
  for (name_id in identifier_count$Var1){
    
    subset_table <- subset(dataframe_calib, identifier == name_id)
    
    subset_bin <- bin_distance_values_new(subset_table$intensityDistanceCalibrated, 
                                          subset_table[[bin_param_calib]],
                                          "binned_calibrated", 
                                          bin_width_calib, 
                                          upper_limit_calib)
    
    
    subset_list[[name_id]] <- subset_bin
    
  }
  
  binned_list <- do.call("rbind", subset_list)
  binned_list1 <- tibble::rownames_to_column(binned_list, "nameindex")
  binned_list1_indices <- str_split_fixed(binned_list1$nameindex, "\\.", 2)
  binned_list2 <- cbind(binned_list1_indices, binned_list1)
  colnames(binned_list2)[1] <- "identifier"
  colnames(binned_list2)[2] <- "index"
  binned_list2$nameindex <- NULL
  head(binned_list2)
  
  plotlist <- list()
  
  # plots intensity profiles
  map_calib <- ggplot(data=binned_list2,  aes(x=row, y=binned_calibrated, color=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() + 
    scale_x_continuous(breaks=binned_list2$row, labels=binned_list2$bin) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    xlab("Normalized distance from nucleus") +
    ggtitle(sprintf("Intensity map distance normalized \n%s channel", channel_string))
  
  ggsave(plot = map_calib,
         file=paste0(plots, .Platform$file.sep, sprintf("%s_intensityMap",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- map_calib
  
  # ------------------------------------------------------------------------------
  # create binned data 
  identifier_count <- as.data.frame(table(dataframe_norm$identifier))
  
  subset_list_norm <- list()
  
  for (name_id in identifier_count$Var1){
    
    subset_table <- subset(dataframe_norm, identifier == name_id)
    
    subset_bin_norm <- bin_distance_values_new(subset_table$intensityDistanceNormalized, 
                                               subset_table[[bin_param_norm]], 
                                               "binned_normalized", 
                                               bin_width_normalized, 
                                               upper_limit_normalized)
    
    subset_list_norm[[name_id]] <- subset_bin_norm
    
  }
  
  binned_list_norm <- do.call("rbind", subset_list_norm)
  binned_list_norm1 <- tibble::rownames_to_column(binned_list_norm, "nameindex")
  binned_list_norm1_indices <- str_split_fixed(binned_list_norm1$nameindex, "\\.", 2)
  binned_list_norm2 <- cbind(binned_list_norm1_indices, binned_list_norm1)
  colnames(binned_list_norm2)[1] <- "identifier"
  colnames(binned_list_norm2)[2] <- "index"
  binned_list_norm2$nameindex <- NULL
  head(binned_list_norm2)
  
  # ------------------------------------------------------------------------------
  # plot of binned data
  # plots distance normalized intensity profiles
  map_norm <- ggplot(data=binned_list_norm2,  aes(x=row, y=binned_normalized, color=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() + 
    scale_x_continuous(breaks=binned_list_norm2$row, labels=binned_list_norm2$bin) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    xlab("Normalized distance from nucleus") +
    ggtitle(sprintf("Intensity map distance normalized \n%s channel", channel_string))
  
  ggsave(plot = map_norm,
         file=paste0(plots, .Platform$file.sep, sprintf("%s_intensityMap_normDist",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- map_norm
  
  # peak normalisation
  name_count_value <- as.data.frame(table(binned_list_norm2$identifier))
  col_intensity_map <- ncol(binned_list_norm2)
  value_list_peak <- list()
  
  for (name in name_count_value$Var1){
    
    data_per_name <- subset(binned_list_norm2, identifier == name)
    
    max_value_profiles = max(data_per_name$binned_normalized, na.rm = TRUE)
    data_per_name$orga_peak_norm <- sapply(data_per_name$binned_normalized, function(x){x /  max_value_profiles})
    
    if ( col_intensity_map == 6) {
      
      max_value_profiles = max(data_per_name$measure_peak_norm, na.rm = TRUE)
      data_per_name$measure_peak_norm <- sapply(data_per_name$measure_peak_norm, function(x){x /  max_value_profiles})
      
    }
    
    value_list_peak[[name]] <- data_per_name
    
  }
  
  # binds collection of normalized value plots and binds them into one dataframe
  norm_list_value <- do.call("rbind", value_list_peak)
  head(norm_list_value)
  
  # plots peak normalized and distance normalized intensity profiles
  map_peak <- ggplot(data=norm_list_value,  aes(x=row, y=orga_peak_norm, color=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() + 
    scale_x_continuous(breaks=norm_list_value$row, labels=norm_list_value$bin) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    xlab("Normalized distance from nucleus") +
    ggtitle(sprintf("Intensity map distance normalized \n%s channel", channel_string))
  
  ggsave(plot = map_peak,
         file=paste0(plots, .Platform$file.sep, sprintf("%s_intensityMap_normDistPeak",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- map_peak
  
  return(plotlist)
  
}