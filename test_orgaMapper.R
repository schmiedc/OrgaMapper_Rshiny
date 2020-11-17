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
source("plot_intensity_ratio.R")
# ==============================================================================
# Params
# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test/"
#directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test_4thChannel/"

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
                                               cal_distance_nucleus,
                                               norm_distance_nucleus,
                                               plot_background_subtract)

do.call(grid.arrange, cell_plots)
do.call(grid.arrange, detection_plots)

if (analyze_signal_profiles) {
  
  # ------------------------------------------------------------------------------
  # collect individual files
  print("Computing individual intensity maps")
  individual_intensity_maps <- collect_individual_profiles_new(directory, 
                                                               series_regex, 
                                                               single_series, 
                                                               cell_measure_filter)
  rownames(individual_intensity_maps) <- c()
  head(individual_intensity_maps)
  
  # create intensity ratio data and plots
  print("Computing and plotting intensity ratio")
  intensity_ratio_results <- compute_intensity_ration(individual_intensity_maps, 
                                                      10, 
                                                      bin_width, 
                                                      0)
  
  plot_intensity_ration(intensity_ratio_results, "orga", plots_intensity)
  
  # group intensity maps
  print("Computing mean of individual intensity maps")
  value_lists <- grouped_intensity_map(individual_intensity_maps)
  
  intensity_map_result <- value_lists$raw
  intensity_map_result_norm <- value_lists$norm
  
  head(intensity_map_result)
  # ------------------------------------------------------------------------------
  print("Saving raw intensity maps")
  write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
             intensity_map_result, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  # ------------------------------------------------------------------------------
  print("Plotting intensity maps")
  
  orga_plots <- plot_intensity_map(intensity_map_result, 
                     intensity_map_result_norm, 
                     "orga", 
                     bin_width, 
                     upper_limit,
                     bin_width_norm,
                     upper_limit_norm,
                     plots_intensity)
  
  do.call(grid.arrange, orga_plots)
  
  if (cell_column == 10 && orga_column == 10) {
    
    measure_plots <- plot_intensity_map(intensity_map_result, 
                                intensity_map_result_norm, 
                                "measure", 
                                bin_width, 
                                upper_limit,
                                bin_width_norm,
                                upper_limit_norm,
                                plots_intensity)
    
    do.call(grid.arrange, measure_plots)
    
  }
    
}

