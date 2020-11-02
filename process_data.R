library(tidyverse)

read_collected_files <- function(inputdir, 
                                 name, 
                                 series_switch, 
                                 series_regex) {
  
  file_name <- paste0(inputdir,name)
  
  file <- read.csv(file_name, header = TRUE)
  
  # TODO needs to default to sensible value if nothing found
  if (series_switch) {
    
    file$Series <- str_extract(file$name, series_regex)
    
    file$Name <- str_remove(file$name, series_regex)

    file$Name <- str_remove(file$name, "(_|-| )($)")

  }
  
  return (file)
  
}

process_cell_measurements <- function(data, 
                                      lower, 
                                      upper, 
                                      column_cell_table, 
                                      column_orga_table,
                                      filter) {
  
  if (filter) {
    
    data_filter <- subset(data, 
                          ferets >= lower & ferets <= upper)
    
  } else {
    
    data_filter <- data
    
  }
  
  # background subtraction for mean organelle intensity per cell
  data_filter$orgaMeanIntensityBacksub <- 
    data_filter$orgaMeanIntensity - data_filter$orgaMeanBackground
  
  if (column_cell_table == 10 && column_orga_table == 10) {
    
    data_filter$measureMeanIntensityBacksub <- 
      data_filter$measureMeanIntensity - data_filter$measureMeanBackground
    
  }
  
  return (data_filter)
  
}

process_orga_measurements <- function(cell_data,
                                      orga_data,
                                      column_cell_table,
                                      column_orga_table) {
  
  merge <- merge(cell_data,
                 orga_data,
                 by = c("name", "series", "cell"))
  
  # background subtraction for detection intensity
  merge$orgaDetectionPeakBacksub <- merge$orgaDetectionPeak - merge$orgaMeanBackground
  
  if (column_cell_table == 10 && column_orga_table == 10) {
    
    merge$measureDetectionPeakBacksub <- merge$measureDetectionPeak - merge$measureMeanBackground
    
  }
  
  merge$detectionDistanceNormalized <- merge$detectionDistanceCalibrated / merge$ferets
  
  return (merge)
  
}

create_summary_table <- function(full_table,
                                 cell_data) {
  
  summary <- full_table %>% 
    group_by(name, series, cell) %>% 
    summarise(across(detectionDistanceRaw:detectionDistanceNormalized, mean, na.rm =TRUE ), .groups = 'drop') %>% 
    rename_at(vars(-name, -series, -cell),function(x) paste0(x,".mean"))
  
  merge <- merge(cell_data,
                 summary, 
                 by = c("name", "series", "cell"))
  
  return (merge)
  
}
