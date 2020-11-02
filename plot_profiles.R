library(ggplot2)
library(tidyverse)
library(lazyeval)

plot_profiles <- function(value_collected,
                          value_norm_collected,
                          channel_string,
                          plots,
                          background_subtract) {
  
  summary_value_norm <- data.frame(Date=as.Date(character()),
                                   File=character(), 
                                   User=character(), 
                                   stringsAsFactors=FALSE) 
  
  summary_value <- data.frame(Date=as.Date(character()),
                              File=character(), 
                              User=character(), 
                              stringsAsFactors=FALSE) 
  
  backsub_norm = ""
  backsub = ""
  raw_norm = ""
  raw = ""
  
  if (channel_string == "organelle") {
    
    backsub_norm = "orgaIntensityBacksub_BinNorm"
    backsub = "orgaIntensityBacksub_Bin"
    raw_norm = "orgaIntensity_BinNorm"
    raw = "orgaIntensity_Bin"
    
  } else if (channel_string == "measure") {
    
    backsub_norm = "measureIntensityBacksub_BinNorm"
    backsub = "measureIntensityBacksub_Bin"
    raw_norm = "measureIntensity_BinNorm"
    raw = "measureIntensity_Bin"
    
  }

  # creates summary data for plots
  if (background_subtract) {
    
    summary_value_norm <- value_norm_collected %>% 
      group_by(identifier, bin) %>% 
      summarise(across(interp(backsub_norm), mean, na.rm =TRUE ))
    
    summary_value <- value_collected %>% 
      group_by(identifier, bin, row) %>% 
      summarise(across(interp(backsub), mean, na.rm =TRUE ))
    
  } else {
    
    summary_value_norm <- value_norm_collected %>% 
      group_by(identifier, bin) %>% 
      summarise(across(interp(raw_norm), mean, na.rm =TRUE ))
    
    summary_value <- value_collected %>% 
      group_by(identifier, bin, row) %>% 
      summarise(across(interp(raw), mean, na.rm =TRUE ))
    
  }
  
  names(summary_value)[4] <- "value"
  names(summary_value_norm)[3] <- "value"

  plotlist <- list()
  
  # plot raw data
  raw_profile <- ggplot(data=summary_value, aes(x=reorder(bin,row), y=value, group=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    ggtitle(sprintf("Intensity profile %s channel", channel_string)) +
    xlab("Distance from nucleus (µm)") 
  
  ggsave(plot = raw_profile,
         file=paste0(plots, .Platform$file.sep, sprintf("profile_%s",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- raw_profile

  profile_norm <- ggplot(data=summary_value_norm, aes(x=bin, y=value, group=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() + 
    aes(x = fct_inorder(bin)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    xlab("Normalized distance from nucleus") +
    ggtitle(sprintf("Norm profile %s channel", channel_string))
  
  ggsave(plot = profile_norm,
         file=paste0(plots, .Platform$file.sep, sprintf("NormProfile_%s",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- profile_norm
  
  # peak normalisation
  name_count_value <- as.data.frame(table(summary_value_norm$identifier))
  value_list <- list()
  
  for (name in name_count_value$Var1){
    
    data_per_name <- subset(summary_value_norm, identifier == name)
    max_value_profiles = max(data_per_name$value, na.rm = TRUE)
    data_per_name$peak_norm <- sapply(data_per_name$value, function(x){x /  max_value_profiles})
    
    value_list[[name]] <- data_per_name
    
  }
  
  # binds collection of normalized value plots and binds them into one dataframe
  norm_list_value <- do.call("rbind", value_list)
  
  # plots peak normalized and distance normalized intensity profiles
  profile_peak <- ggplot(data=norm_list_value, aes(x=bin, y=peak_norm, group=identifier)) +
    geom_line(aes(color=identifier), na.rm=TRUE) +
    geom_point(aes(color=identifier), na.rm=TRUE) +
    lineplot_theme() + 
    aes(x = fct_inorder(bin)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Normalized fluorescent intensity") +
    xlab("Normalized distance from Nucleus") +
    ggtitle(sprintf("Peak norm profile %s channel", channel_string))
  
  ggsave(plot = profile_peak,
         file=paste0(plots, .Platform$file.sep, sprintf("PeakNormProfile_%s",channel_string), ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  plotlist[[length(plotlist) + 1 ]] <- profile_peak
  
  return (plotlist)
}