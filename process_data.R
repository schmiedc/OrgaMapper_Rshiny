library(tidyverse)

read_collected_files <- function(inputdir, 
                                 name, 
                                 series_switch, 
                                 series_regex) {
  
  file_name <- paste0(inputdir,name)
  
  file <- read.csv(file_name, header = TRUE)
  
  # TODO needs to default to sensible value if nothing found
  if (series_switch) {
    
    file$Series <- str_extract(file$Name, series_regex)
    
    file$Name <- str_remove(file$Name, series_regex)

    file$Name <- str_remove(file$Name, "(_|-| )($)")

  }
  
  return (file)
  
}

process_cell_measurements <- function(data, 
                                      lower, 
                                      upper, 
                                      column_cell_table, 
                                      column_orga_table) {
  
  data_filter <- subset(data, 
                        Ferets >= lower & Ferets <= upper)
  
  # background subtraction for mean organelle intensity per cell
  data_filter$MeanOrgaBackSub <- 
    data_filter$MeanValueOrga - data_filter$MeanBackgroundOrga
  
  if (column_cell_table == 10 && column_orga_table == 10) {
    
    data_filter$MeanMeasureBackSub <- 
      data_filter$MeanValueMeasure - data_filter$MeanBackgroundMeasure
    
  }
  
  return (data_filter)
  
}

process_orga_measurements <- function(cell_data,
                                      orga_data,
                                      column_cell_table,
                                      column_orga_table) {
  
  merge <- merge(cell_data,
                 orga_data,
                 by = c("Name", "Series", "Cell"))
  
  # background subtraction for detection intensity
  merge$PeakDetectBackSub <- merge$PeakDetectionInt - merge$MeanBackgroundOrga
  
  if (column_cell_table == 10 && column_orga_table == 10) {
    
    merge$PeakMeasureBackSub <- merge$PeakMeasureInt - merge$MeanBackgroundMeasure
    
  }
  
  merge$DistanceNorm <- merge$DistanceCal / merge$Ferets
  
  return (merge)
  
}

create_summary_table <- function(full_table,
                                 cell_data) {
  
  summary <- full_table %>% 
    group_by(Name, Series, Cell) %>% 
    summarise(across(DistanceRaw:DistanceNorm, mean, na.rm =TRUE ), .groups = 'drop') %>% 
    rename_at(vars(-Name, -Series, -Cell),function(x) paste0(x,".mean"))
  
  merge <- merge(cell_data,
                 summary, 
                 by = c("Name", "Series", "Cell"))
  
  return (merge)
  
}
