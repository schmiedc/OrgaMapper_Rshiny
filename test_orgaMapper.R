setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")

packages <- c("shiny", "shinyFiles", "openxlsx", "ggplot2", "gridExtra", "tidyverse", "lazyeval")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library("openxlsx")
library(gridExtra)
library(lazyeval)
source("process_data.R")
source("plot_data.R")
source("process_profiles.R")
source("plot_profiles.R")
source("plot_intensity_ratio.R")
# ==============================================================================
# Params
# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test/"
# directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test_4thChannel/"

result_name = "Analysis_test"

# filter for feret's diameter
feret_filter = TRUE
feret_lower = 0
feret_upper = 600

# determine range for plots
cal_distance_nucleus = 75
norm_distance_nucleus = 0.7

# TODO if file contains series number or the already present column
# needs to default to something sensible if not possible
single_series = FALSE
series_regex = "(?<=_)\\d*($)"

# TODO apply background subtraction for plots
plot_background_subtract = TRUE

# analyze signal profiles
analyze_signal_profiles = FALSE

# Binning for intensity profiles
# or different method for binning
upper_limit_norm = 1
bin_width_norm = 0.05

upper_limit = 75
bin_width = 2

# ==============================================================================
# where to save the data
out_dir =  directory
result_path <- file.path(out_dir, result_name, fsep = .Platform$file.sep)

# plot dir
plots_distance <- file.path(out_dir, "plot_distance_map", fsep = .Platform$file.sep)
dir.create(plots_distance, showWarnings = FALSE)

# create directory for intensity maps
plots_intensity <- file.path(out_dir, "plot_intensity_map", fsep = .Platform$file.sep)
dir.create(plots_intensity, showWarnings = FALSE)

# ==============================================================================
name_distance = "organelleDistance.csv"
name_cell_measure = "cellMeasurements.csv"

organelle_distance <- read_collected_files(directory, 
                                           name_distance, 
                                           single_series, 
                                           series_regex)

cell_measure <- read_collected_files(directory, 
                                     name_cell_measure, 
                                     single_series, 
                                     series_regex)

cell_column <- ncol(cell_measure)
orga_column <- ncol(organelle_distance)


cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                 feret_lower, 
                                                 feret_upper,
                                                 cell_column,
                                                 orga_column,
                                                 feret_filter)

merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                  organelle_distance,
                                                  cell_column,
                                                  orga_column)

merged_summary <- create_summary_table(merge_cell_organelle,
                                       cell_measure_filter)

# ------------------------------------------------------------------------------
# renaming for organelle result tables
if (cell_column == 10 && orga_column == 10) {
  
  merge_cell_organelle_result <- merge_cell_organelle %>%
    rename(
      cell_area = cellArea,
      numberOfDetections = numberDetections,
      orga_intensity = orgaMeanIntensity,
      orga_background = orgaMeanBackground,
      measure_intensity = measureMeanIntensity,
      measure_background = measureMeanBackground,
      orga_intensity_backsub = orgaMeanIntensityBacksub,
      measure_intensity_backsub = measureMeanIntensityBacksub,
      x_detection = xDetection,
      y_detection = yDetection,
      orga_distance_pixel = detectionDistanceRaw,
      orga_distance_calibrated = detectionDistanceCalibrated,
      orga_detection_peak = orgaDetectionPeak,
      measure_detection_peak = measureDetectionPeak,
      orga_detection_peak_backsub = orgaDetectionPeakBacksub,
      measure_detection_peak_backsub = measureDetectionPeakBacksub,
      orga_distance_normalized = detectionDistanceNormalized
           )
  
} else {
  
  merge_cell_organelle_result <- merge_cell_organelle %>%
    rename(
      cell_area = cellArea,
      numberOfDetections = numberDetections,
      orga_intensity = orgaMeanIntensity,
      orga_background = orgaMeanBackground,
      orga_intensity_backsub = orgaMeanIntensityBacksub,
      x_detection = xDetection,
      y_detection = yDetection,
      orga_distance_pixel = detectionDistanceRaw,
      orga_distance_calibrated = detectionDistanceCalibrated,
      orga_detection_peak = orgaDetectionPeak,
      orga_detection_peak_backsub = orgaDetectionPeakBacksub,
      orga_distance_normalized = detectionDistanceNormalized
    )
  
}

# ------------------------------------------------------------------------------
# save processed data
write.xlsx(file = paste0( result_path, "_detection.xlsx", sep = ""), 
           merge_cell_organelle_result, 
           sheetName="Sheet1",  
           col.names=TRUE, 
           row.names=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# ------------------------------------------------------------------------------

if (cell_column == 10 && orga_column == 10) {
  
  merged_summary_result <- merged_summary %>%
    rename(
      cell_area = cellArea,
      orga_numberOfDetections = numberDetections,
      orga_intensity = orgaMeanIntensity,
      orga_background = orgaMeanBackground,
      measure_intensity = measureMeanIntensity,
      measure_background = measureMeanBackground,
      orga_intensity_backsub = orgaMeanIntensityBacksub,
      measure_intensity_backsub = measureMeanIntensityBacksub,
      orga_meanDistance_pixel = detectionDistanceRaw.mean,
      orga_meanDistance_calibrated = detectionDistanceCalibrated.mean,
      orga_intensityOnDetection = orgaDetectionPeak.mean,
      measure_intensityOnDetection = measureDetectionPeak.mean,
      orga_intensityOnDetection_backsub = orgaDetectionPeakBacksub.mean,
      measure_intensityOnDetection_backsub = measureDetectionPeakBacksub.mean,
      orga_meanDistance_normalized = detectionDistanceNormalized.mean
      )
  
} else {
  
  merged_summary_result <- merged_summary %>%
    rename(
      cell_area = cellArea,
      orga_numberOfDetections = numberDetections,
      orga_intensity = orgaMeanIntensity,
      orga_background = orgaMeanBackground,
      orga_intensity_backsub = orgaMeanIntensityBacksub,
      orga_meanDistance_pixel = detectionDistanceRaw.mean,
      orga_meanDistance_calibrated = detectionDistanceCalibrated.mean,
      orga_intensityOnDetection = orgaDetectionPeak.mean,
      orga_intensityOnDetection_backsub = orgaDetectionPeakBacksub.mean,
      orga_meanDistance_normalized = detectionDistanceNormalized.mean
    )
  
}
# merged_summary
write.xlsx(file = paste0( result_path,  "_cell.xlsx", sep = ""), 
           merged_summary_result, 
           sheetName="Sheet1",  
           col.names=TRUE, 
           row.names=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# ------------------------------------------------------------------------------
# plot data
cell_plots <- plot_cell_measurements(cell_measure_filter,
                                     plots_distance,
                                     cell_column,
                                     orga_column,
                                     plot_background_subtract)



detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                               merged_summary,
                                               plots_distance,
                                               cell_column,
                                               orga_column,
                                               cal_distance_nucleus,
                                               norm_distance_nucleus,
                                               plot_background_subtract)

do.call(grid.arrange, cell_plots)
do.call(grid.arrange, detection_plots)

# ==============================================================================

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

# ------------------------------------------------------------------------------
# collect individual files
individual_intensity_maps <- collect_individual_profiles_new(directory, series_regex, single_series, cell_measure_filter)
rownames(individual_intensity_maps) <- c()
head(individual_intensity_maps)

# create intensity ratio data and plots
intensity_ratio_results <- compute_intensity_ration(individual_intensity_maps, 10, bin_width, 0)
plot_intensity_ration(intensity_ratio_results, "orga", plots_intensity)

# group intensity maps
value_lists <- grouped_intensity_map(individual_intensity_maps)

head(value_lists$raw)
head(value_lists$norm)


# ------------------------------------------------------------------------------
# create binned data 
identifier_count <- as.data.frame(table(value_list_treat_norm$identifier))

subset_list <- list()

for (name_id in identifier_count$Var1){

  subset_table <- subset(value_list_treat_norm, identifier == name_id)
  
  subset_bin <- bin_distance_values_new(subset_table$intensityDistanceNormalized, 
                                        subset_table$mean, 
                                        "binned_values", 
                                        bin_width_norm, 
                                        upper_limit_norm)
  
  subset_list[[name_id]] <- subset_bin
  
}

binned_value_list <- do.call("rbind", subset_list)
binned_value_list1 <- tibble::rownames_to_column(binned_value_list, "nameindex")
binned_value_list1_indices <- str_split_fixed(binned_value_list1$nameindex, "\\.", 2)
binned_value_list2 <- cbind(binned_value_list1_indices, binned_value_list1)
colnames(binned_value_list2)[1] <- "identifier"
colnames(binned_value_list2)[2] <- "index"
head(binned_value_list2)

# plots peak normalized and distance normalized intensity profiles
ggplot(data=binned_value_list2,  aes(x=row, y=binned_values, color=identifier)) +
  geom_line(aes(color=identifier), na.rm=TRUE) +
  geom_point(aes(color=identifier), na.rm=TRUE) +
  ylab("Fluorescent intensity (A.U.)") +
  xlab("Normalized distance from nucleus") +
  ggtitle(sprintf("Intensity map distance normalized \nOrga channel"))


















# ==============================================================================
if (analyze_signal_profiles) {
  
  name_value_measure = "intensityDistance.csv"
  
  profile_collected <- collect_individual_profiles(directory,
                                                   name_value_measure,
                                                   cell_measure_filter,
                                                   single_series,
                                                   series_regex,
                                                   upper_limit,
                                                   bin_width,
                                                   upper_limit_norm,
                                                   bin_width_norm)
  
  value_list <- profile_collected$raw
  value_list_norm <- profile_collected$norm
  rownames(value_list) <- c()
  rownames(value_list_norm) <- c()
  
  head(value_list)
  head(value_list_norm)
  # ------------------------------------------------------------------------------
  if (cell_column == 10 && orga_column == 10) {
  
    value_list_result <- value_list %>%
      rename(
        cell_area = cellArea,
        orga_numberOfDetections = numberDetections,
        orga_intensity = orgaMeanIntensity,
        orga_background = orgaMeanBackground,
        measure_intensity = measureMeanIntensity,
        measure_background = measureMeanBackground,
        orga_intensity_backsub = orgaMeanIntensityBacksub,
        measure_intensity_backsub = measureMeanIntensityBacksub,
        orga_intensityMap_backsub = orgaIntensityBacksub_Bin,
        orga_intensityMap = orgaIntensity_Bin,
        measure_intensityMap_backsub = measureIntensityBacksub_Bin,
        measure_intensityMap = measureIntensity_Bin
      )
    
    } else {
      
      value_list_result <- value_list %>%
        rename(
          cell_area = cellArea,
          orga_numberOfDetections = numberDetections,
          orga_intensity = orgaMeanIntensity,
          orga_background = orgaMeanBackground,
          orga_intensity_backsub = orgaMeanIntensityBacksub,
          orga_intensityMap_backsub = orgaIntensityBacksub_Bin,
          orga_intensityMap = orgaIntensity_Bin
        )
      
    }
  
  write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
             value_list_result, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  # ----------------------------------------------------------------------------
  # compute & plot intensity ratio
  intensity_ratio_results <- compute_intensity_ration(value_list_result, 10, bin_width, 0)

  plot_intensity_ration(intensity_ratio_results, "orga", plots_intensity)
  
  # ----------------------------------------------------------------------------
  if (cell_column == 10 && orga_column == 10) {
    
    value_list_norm_result <- value_list_norm %>%
      rename(
        cell_area = cellArea,
        orga_numberOfDetections = numberDetections,
        orga_intensity = orgaMeanIntensity,
        orga_background = orgaMeanBackground,
        measure_intensity = measureMeanIntensity,
        measure_background = measureMeanBackground,
        orga_intensity_backsub = orgaMeanIntensityBacksub,
        measure_intensity_backsub = measureMeanIntensityBacksub,
        orga_intensityMapNorm_backsub = orgaIntensityBacksub_BinNorm,
        orga_intensityMapNorm = orgaIntensity_BinNorm,
        measure_intensityMapNorm_backsub = measureIntensityBacksub_BinNorm,
        measure_intensityMapNorm = measureIntensity_BinNorm
      )
    
  } else {
    
    value_list_norm_result <- value_list_norm %>%
      rename(
        cell_area = cellArea,
        orga_numberOfDetections = numberDetections,
        orga_intensity = orgaMeanIntensity,
        orga_background = orgaMeanBackground,
        orga_intensity_backsub = orgaMeanIntensityBacksub,
        orga_intensityMapNorm_backsub = orgaIntensityBacksub_BinNorm,
        orga_intensityMapNorm = orgaIntensity_BinNorm
        
      )
    
  }
  
  write.xlsx(file = paste0( result_path,  "_intensityProfile_norm.xlsx", sep = ""), 
             value_list_norm_result, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)

  # ----------------------------------------------------------------------------
  organelle_profile <- plot_profiles(value_list, 
                                value_list_norm, 
                                "orga", 
                                plots_intensity, 
                                plot_background_subtract)
  
  if (cell_column == 10 && orga_column == 10) {
    
    measure_profile <- plot_profiles(value_list, 
                                      value_list_norm, 
                                      "measure", 
                                      plots_intensity, 
                                      plot_background_subtract)
    
    do.call(grid.arrange, measure_profile)
    
  }
  
  do.call(grid.arrange, organelle_profile)
  
}

