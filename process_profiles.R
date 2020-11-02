library(tidyverse)

# creates binned data
bin_distance_values <- function(bin, value, variable_name, width, limit) {
  
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

process_profile_data <- function(value_data,
                                 per_cell_data,
                                 c_col,
                                 v_col) {
  
  # merge cell measurements with intensity profiles
  merge_table <- merge(per_cell_data, 
                             value_data, 
                             by = c("identifier", "series", "cell"))
  
  head(merge_table)
  # distance normalization
  merge_table$intensityDistanceNormalized <- merge_table$intensityDistanceCalibrated / merge_table$ferets
  
  # background subtraction for detection intensity
  merge_table$orgaIntensityBacksub <- merge_table$orgaIntensity - merge_table$orgaMeanBackground
  
  if (c_col == 12 && v_col == 9) {
    
    merge_table$measureIntensityBackSub <- merge_table$measureIntensity - merge_table$measureMeanBackground
    
  }
  
  return (merge_table)
  
}

collect_individual_profiles <- function(inputdir,
                                        name,
                                        cell_measure_data,
                                        series,
                                        regular_expression,
                                        limit,
                                        width,
                                        limit_norm,
                                        width_norm) {

  cell_col = ncol(cell_measure_data)

  file_name <- paste0(inputdir,name)
  
  value_filenames <- list.files(path = inputdir,
                                recursive=TRUE, 
                                pattern=name,
                                full.names = TRUE)

  value_list_norm <- list()
  value_list <- list()
  
  print("Retrieving individual files")
  
  for (name in value_filenames) {
    
    print(paste("Processing file: ", name))
    
    value_measure  <- read.csv(name, header = TRUE)
    
    # check if intDistance.csv is empty if true skip
    if(nrow(value_measure) == 0) {
      print(paste("Skipping file: ", name))
      next
    }
    
    # TODO needs to default to something that is sensible if invalid
    if (series) {
      
      value_measure$series <- str_extract(value_measure$identifier, regular_expression)
      
      value_measure$identifier <- str_remove(value_measure$identifier, regular_expression)
      
      # removes trailing underscore or hypen
      value_measure$identifier <- str_remove(value_measure$identifier, "(_|-| )($)")
      
    }
    
    value_col = ncol(value_measure)
    
    print("Merging with cell data")
    merge_table_value <- process_profile_data(value_measure,
                                              cell_measure_data,
                                              cell_col,
                                              value_col)

    print("Bin normalized distance of intensity background subtracted values")
    value_result_norm <- bin_distance_values(merge_table_value$intensityDistanceNormalized, 
                                             merge_table_value$orgaIntensityBacksub, 
                                             "backsub_bin_orga_norm",
                                             width_norm,
                                             limit_norm)
    
    print("Bin raw distance of intensity background subtracted values")
    value_result <- bin_distance_values(merge_table_value$intensityDistanceCalibrated, 
                                        merge_table_value$orgaIntensityBacksub, 
                                        "backsub_bin_orga",
                                        width,
                                        limit)
    
    print("Bin normalized distance of intensity raw values")
    value_result_norm_raw <- bin_distance_values(merge_table_value$intensityDistanceNormalized, 
                                                 merge_table_value$orgaIntensity, 
                                                 "raw_bin_orga_norm",
                                                 width_norm,
                                                 limit_norm)
    
    print("Bin raw distance of intensity raw values")
    value_result_raw <- bin_distance_values(merge_table_value$intensityDistanceCalibrated, 
                                            merge_table_value$orgaIntensity, 
                                            "raw_bin_orga",
                                            width,
                                            limit)
    
    value_result_norm <- merge(value_result_norm, value_result_norm_raw, by = c("bin","row"))
    value_result <- merge(value_result, value_result_raw, by = c("bin","row"))
    
    if (cell_col == 12 && value_col == 9) {
      
      print("Mapping measurement channel")
      binned_measure_value_norm <- bin_distance_values(merge_table_value$intensityDistanceNormalized, 
                                                       merge_table_value$measureIntensityBackSub, 
                                                       "backsub_bin_measure_norm",
                                                       width_norm,
                                                       limit_norm)
      
      binned_measure_value <- bin_distance_values(merge_table_value$intensityDistanceCalibrated, 
                                                  merge_table_value$measureIntensityBackSub, 
                                                  "backsub_bin_measure",
                                                  width,
                                                  limit)
      
      binned_measure_value_norm_raw <- bin_distance_values(merge_table_value$intensityDistanceNormalized, 
                                                           merge_table_value$measureIntensity, 
                                                           "raw_bin_measure_norm",
                                                           width_norm,
                                                           limit_norm)
      
      binned_measure_value_raw <- bin_distance_values(merge_table_value$intensityDistanceCalibrated, 
                                                      merge_table_value$measureIntensity, 
                                                      "raw_bin_measure",
                                                      width,
                                                      limit)
      
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