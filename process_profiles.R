library(tidyverse)

# creates binned data
bin_distance_values <- function(bin, value, variable_name, bin_width, upper_limit) {
  
  output_apply <- tapply(value, 
                         cut(bin, seq(0, upper_limit, by=bin_width)), 
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

process_profile_data <- function(value_data,
                                 per_cell_data,
                                 c_col,
                                 v_col) {
  
  
  
  # merge cell measurements with intensity profiles
  merge_table <- merge(per_cell_data, 
                             value_data, 
                             by = c("Name", "Series", "Cell"))
  
  # distance normalization
  merge_table$DistanceNorm <- merge_table$DistanceCal / merge_table$Ferets
  
  # background subtraction for detection intensity
  merge_table$orgaIntBackSub <- merge_table$orgaInt - merge_table$MeanBackgroundOrga
  
  if (c_col == 10 && v_col == 9) {
    
    merge_table$measureIntBackSub <- merge_table$measureInt - merge_table$MeanBackgroundMeasure
    
  }
  
  return (merge_table)
  
}

collect_individual_profiles <- function(inputdir,
                                        name,
                                        cell_measure_data,
                                        series,
                                        regular_expression) {

  cell_col = ncol(cell_measure_data)
  
  file_name <- paste0(inputdir,name)
  
  value_filenames <- list.files(path = inputdir,
                                recursive=TRUE, 
                                pattern=name,
                                full.names = TRUE)

  value_list_norm <- list()
  value_list <- list()
  
  for (name in value_filenames) {
    
    value_measure  <- read.csv(name, header = TRUE)
    
    # TODO needs to default to something that is sensible if invalid
    if (series) {
      
      value_measure$Series <- str_extract(value_measure$Name, series_regex)
      
      value_measure$Name <- str_remove(value_measure$Name, series_regex)
      
      # removes trailing underscore or hypen
      value_measure$Name <- str_remove(value_measure$Name, regular_expression)
      
    }
    
    value_col = ncol(value_measure)
    
    merge_table_value <- process_profile_data(value_measure,
                                              cell_measure_data,
                                              cell_col,
                                              value_col)
    
    value_result_norm <- bin_distance_values(merge_table_value$DistanceNorm, 
                                             merge_table_value$orgaIntBackSub, 
                                             "backsub_bin_orga_norm",
                                             bin_width_norm,
                                             upper_limit_norm)
    
    value_result <- bin_distance_values(merge_table_value$DistanceCal, 
                                        merge_table_value$orgaIntBackSub, 
                                        "backsub_bin_orga",
                                        bin_width,
                                        upper_limit)
    
    value_result_norm_raw <- bin_distance_values(merge_table_value$DistanceNorm, 
                                                 merge_table_value$orgaInt, 
                                                 "raw_bin_orga_norm",
                                                 bin_width_norm,
                                                 upper_limit_norm)
    
    value_result_raw <- bin_distance_values(merge_table_value$DistanceCal, 
                                            merge_table_value$orgaInt, 
                                            "raw_bin_orga",
                                            bin_width,
                                            upper_limit)
    
    value_result_norm <- merge(value_result_norm, value_result_norm_raw, by = c("bin","row"))
    value_result <- merge(value_result, value_result_raw, by = c("bin","row"))
    
    if (cell_col == 10 && value_col == 9) {
      
      binned_measure_value_norm <- bin_distance_values(merge_table_value$DistanceNorm, 
                                                       merge_table_value$measureIntBackSub, 
                                                       "backsub_bin_measure_norm",
                                                       bin_width_norm,
                                                       upper_limit_norm)
      
      binned_measure_value <- bin_distance_values(merge_table_value$DistanceCal, 
                                                  merge_table_value$measureIntBackSub, 
                                                  "backsub_bin_measure",
                                                  bin_width,
                                                  upper_limit)
      
      binned_measure_value_norm_raw <- bin_distance_values(merge_table_value$DistanceNorm, 
                                                           merge_table_value$measureInt, 
                                                           "raw_bin_measure_norm",
                                                           bin_width_norm,
                                                           upper_limit_norm)
      
      binned_measure_value_raw <- bin_distance_values(merge_table_value$DistanceCal, 
                                                      merge_table_value$measureInt, 
                                                      "raw_bin_measure",
                                                      bin_width,
                                                      upper_limit)
      
      value_result_norm <- merge(value_result_norm, binned_measure_value_norm, by = c("bin","row"))
      value_result <- merge(value_result, binned_measure_value, by = c("bin","row"))
      value_result_norm <- merge(value_result_norm, binned_measure_value_norm_raw, by = c("bin","row"))
      value_result <- merge(value_result, binned_measure_value_raw, by = c("bin","row"))
      
    }
    
    table.info <- merge_table_value[1,1:12]
    value_result_norm <- merge(table.info, value_result_norm, by=NULL)
    value_result<- merge(table.info, value_result, by=NULL)
    
    value_list_norm[[name]] <- value_result_norm
    value_list[[name]] <- value_result
    
  }
  
  # binds collection 
  value_list_norm_coll <- do.call("rbind", value_list_norm)
  value_list_coll <- do.call("rbind", value_list)

  return (list("raw" = value_list_coll, "norm" = value_list_norm_coll ))
  
}