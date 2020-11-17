
# opening files that contain intensity information
collect_individual_profiles_new <- function(inputdir, regular_expression, series, cell_measure_data) {
  
  cell_col = ncol(cell_measure_data)
  
  
  
  name = "intensityDistance.csv"
  
  print("Retrieving individual files")
  
  file_name <- paste0(inputdir,name)
  
  value_filenames <- list.files(path = inputdir,
                                recursive=TRUE, 
                                pattern=name,
                                full.names = TRUE)
  
  value_list <- list()
  
  for (name in value_filenames) {
    
    print(paste("Processing file: ", name))
    
    value_measure  <- read.csv(name, header = TRUE)
    
    value_col = ncol(value_measure)
    
    # check if intDistance.csv is empty if true skip
    if(nrow(value_measure) == 0) {
      print(paste("Skipping file: ", name))
      next
    }
    
    # TODO needs to default to something that is sensible if invalid
    if (series) {
      
      value_measure$series <- str_extract(value_measure$identifier, regular_expression)
      
      value_measure$identifier <- str_remove(value_measure$identifier, regular_expression)
      
      # removes trailing underscore or hyphen
      value_measure$identifier <- str_remove(value_measure$identifier, "(_|-| )($)")
      
    }
    
    # TODO check if this works
    if (cell_col == 12 && value_col == 9) {
      
      value_measure_mean <- value_measure %>% 
        group_by(identifier, series, cell, intensityDistanceCalibrated) %>% 
        summarise(mean_orgaIntensity = mean(orgaIntensity), mean_measureIntensity = mean(measureIntensity))
      
    } else {
      
      value_measure_mean <- value_measure %>% 
        group_by(identifier, series, cell, intensityDistanceCalibrated) %>% 
        summarise(mean_orgaIntensity = mean(orgaIntensity))
      
    }
    
    # merge cell measurements with intensity profiles
    merge_table <- merge(cell_measure_data, 
                         value_measure_mean, 
                         by = c("identifier", "series", "cell"))
    
    # distance normalization
    merge_table$intensityDistanceNormalized <- merge_table$intensityDistanceCalibrated / merge_table$ferets
    
    # background subtraction for detection intensity
    merge_table$orgaIntensityBacksub <- merge_table$mean_orgaIntensity - merge_table$orgaMeanBackground
    
    if (cell_col == 12 && value_col == 9) {
      
      merge_table$measureIntensityBackSub <- merge_table$mean_measureIntensity - merge_table$measureMeanBackground
      
    }
    
    value_list[[name]] <- merge_table
    
  }
  
  value_list_coll <- do.call("rbind", value_list)
  gc()
  return(value_list_coll)
}

bin_distance_values_new <- function(bin, value, variable_name, width, limit) {
  
  output_apply <- tapply(value, 
                         cut(bin, seq(0, limit, by=width)), 
                         mean)
  
  # transform output of tapply to data frame
  binned_values <- as.data.frame(as.table(output_apply))
  colnames(binned_values) <- c("bin", variable_name)
  
  # cleanup data frame
  binned_values$bin <- str_remove_all(binned_values$bin, "[]\\(]")
  binned_values$bin <- str_replace(binned_values$bin, ",", "-")
  
  binned_values <- binned_values %>% mutate(row = row_number())
  
  return(binned_values)
  
}

# create different grouped means
grouped_intensity_map <- function(individual_maps) {
  
  intensity_maps_col = ncol(individual_maps)
  
  if (intensity_maps_col == 18) {
    
    if (plot_background_subtract) {
      
      value_list_treat <- individual_maps %>% 
        group_by(identifier,intensityDistanceCalibrated) %>% 
        summarise(orga_mean = mean(orgaIntensityBacksub), measure_mean = mean(measureIntensityBackSub))
      
      value_list_treat_norm <- individual_maps %>% 
        group_by(identifier, intensityDistanceNormalized) %>% 
        summarise(orga_mean = mean(orgaIntensityBacksub), measure_mean = mean(measureIntensityBackSub))
      
    } else {
      
      value_list_treat <- individual_maps %>% 
        group_by(identifier,intensityDistanceCalibrated) %>% 
        summarise(orga_mean = mean(mean_orgaIntensity), measure_mean = mean(mean_measureIntensity))
      
      value_list_treat_norm <- individual_maps %>% 
        group_by(identifier, intensityDistanceNormalized) %>% 
        summarise(orga_mean = mean(mean_orgaIntensity), measure_mean = mean(mean_measureIntensity))
    }
    
  } else {
    
    if (plot_background_subtract) {
      
      value_list_treat <- individual_maps %>% 
        group_by(identifier,intensityDistanceCalibrated) %>% 
        summarise(orga_mean = mean(orgaIntensityBacksub))
      
      value_list_treat_norm <- individual_maps %>% 
        group_by(identifier, intensityDistanceNormalized) %>% 
        summarise(orga_mean = mean(orgaIntensityBacksub))
      
    } else {
      
      value_list_treat <- individual_maps %>% 
        group_by(identifier,intensityDistanceCalibrated) %>% 
        summarise(orga_mean = mean(mean_orgaIntensity))
      
      value_list_treat_norm <- individual_maps %>% 
        group_by(identifier, intensityDistanceNormalized) %>% 
        summarise(orga_mean = mean(mean_orgaIntensity))
      
    }
    
  }
  
  # peak normalisation
  name_count_value <- as.data.frame(table(value_list_treat_norm$identifier))
  col_intensity_map <- ncol(value_list_treat_norm)
  value_list_peak <- list()
  
  for (name in name_count_value$Var1){
    
    data_per_name <- subset(value_list_treat_norm, identifier == name)
    
    max_value_profiles = max(data_per_name$orga_mean, na.rm = TRUE)
    data_per_name$orga_peak_norm <- sapply(data_per_name$orga_mean, function(x){x /  max_value_profiles})
    
    if ( col_intensity_map == 4) {
      
      max_value_profiles = max(data_per_name$measure_mean, na.rm = TRUE)
      data_per_name$measure_peak_norm <- sapply(data_per_name$measure_mean, function(x){x /  max_value_profiles})
      
    }
    
    value_list_peak[[name]] <- data_per_name
    
  }
  
  # binds collection of normalized value plots and binds them into one dataframe
  norm_list_value <- do.call("rbind", value_list_peak)
  
  return (list("raw" = value_list_treat, "norm" = norm_list_value ))
  
}
