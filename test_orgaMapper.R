setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")

packages <- c("shiny", "shinyFiles", "openxlsx", "ggplot2", "gridExtra", "tidyverse", "lazyeval")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library("openxlsx")
library(gridExtra)
source("process_data.R")
source("plot_data.R")
source("process_profiles.R")
source("plot_profiles.R")
# ==============================================================================
# Params
# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/Test/test_nd2/output_noMeasure/"

result_name = "Analysis_test"

# filter for feret's diameter
feret_filter = TRUE
feret_lower = 0
feret_upper = 600

# determine range for plots
norm_distance_nucleus = 0.7

# TODO if file contains series number or the already present column
# needs to default to something sensible if not possible
single_series = FALSE
series_regex = "(?<=_)\\d*($)"

# TODO apply background subtraction for plots
plot_background_subtract = TRUE

# analyze signal profiles
analyze_signal_profiles = TRUE

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
                                               norm_distance_nucleus,
                                               plot_background_subtract)

do.call(grid.arrange, cell_plots)
do.call(grid.arrange, detection_plots)

# if (analyze_signal_profiles) {
  
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
  
# }

