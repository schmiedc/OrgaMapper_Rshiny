setwd("/home/schmiedc/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")

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
#directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test/"
#directory = "/home/schmiedc/Desktop/OrgaMapper_Data/siArl8b_vs_scr/output_test_4thChannel/"
directory = "/home/schmiedc/Downloads/download/Output/"

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
single_series = TRUE
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

cell_measure <- read_collected_files(directory, 
                                     name_cell_measure, 
                                     single_series, 
                                     series_regex)

organelle_distance <- read_collected_files(directory, 
                                           name_distance, 
                                           single_series, 
                                           series_regex)

# Checks if there were measurements in measurement channel
measureChannelCell = "measureMeanIntensity" %in% colnames(cell_measure)
measureChannelOrganelle = "measureDetectionPeak" %in% colnames(organelle_distance)

cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                 feret_lower, 
                                                 feret_upper,
                                                 measureChannelCell,
                                                 measureChannelOrganelle,
                                                 feret_filter)

merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                  organelle_distance,
                                                  measureChannelCell,
                                                  measureChannelOrganelle)

merged_summary <- create_summary_table(merge_cell_organelle,
                                       cell_measure_filter)

# ------------------------------------------------------------------------------
# saving data

# renaming for organelle result tables
detection_lookup <- c(cell_area = "cellArea",
                      numberOfDetections = "numberDetections",
                      orga_intensity = "orgaMeanIntensity",
                      orga_background = "orgaMeanBackground",
                      x_nucleus_center_mass = "nucleusCenterMassX",
                      y_nucleus_center_mass = "nucleusCenterMassY",
                      measure_intensity = "measureMeanIntensity",
                      measure_background = "measureMeanBackground",
                      orga_intensity_backsub = "orgaMeanIntensityBacksub",
                      measure_intensity_backsub = "measureMeanIntensityBacksub",
                      x_detection = "xDetection",
                      y_detection = "yDetection",
                      orga_distance_nucleus_pixel = "detectionDistanceRaw",
                      orga_distance_nucleus_calibrated = "detectionDistanceCalibrated",
                      orga_detection_peak = "orgaDetectionPeak",
                      measure_detection_peak = "measureDetectionPeak",
                      orga_detection_peak_backsub = "orgaDetectionPeakBacksub",
                      measure_detection_peak_backsub = "measureDetectionPeakBacksub",
                      orga_distance_nucleus_normalized = "detectionDistanceNormalized")

merge_cell_organelle_result <- merge_cell_organelle %>%
  rename(
    any_of(
      
      detection_lookup
      
    )
  )

# TODO: Cells with no detections are now empty
# save processed data
write.xlsx(file = paste0( result_path, "_detection.xlsx", sep = ""), 
           merge_cell_organelle_result, 
           sheetName="Sheet1",  
           colNames=TRUE, 
           rowNames=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# renaming for cell results
cell_lookup <- c(cell_area = "cellArea",
                 orga_numberOfDetections = "numberDetections",
                 orga_intensity = "orgaMeanIntensity",
                 orga_background = "orgaMeanBackground",
                 measure_intensity = "measureMeanIntensity",
                 measure_background = "measureMeanBackground",
                 x_nucleus_center_mass = "nucleusCenterMassX",
                 y_nucleus_center_mass = "nucleusCenterMassY",
                 orga_intensity_backsub = "orgaMeanIntensityBacksub",
                 measure_intensity_backsub = "measureMeanIntensityBacksub",
                 orga_meanDistance_nucleus_pixel = "detectionDistanceRaw.mean",
                 orga_meanDistance_nucleus_calibrated = "detectionDistanceCalibrated.mean",
                 measure_intensityOnDetection = "orgaDetectionPeak.mean",
                 orga_intensityOnDetection_backsub = "orgaDetectionPeakBacksub.mean",
                 measure_intensityOnDetection_backsub = "measureDetectionPeakBacksub.mean",
                 orga_meanDistance_nucleus_normalized = "detectionDistanceNormalized.mean")

merged_summary_result <- merged_summary %>% 
  rename(
    any_of(
      cell_lookup
    )
  )

# merged_summary
write.xlsx(file = paste0( result_path,  "_cell.xlsx", sep = ""), 
           merged_summary_result, 
           sheetName="Sheet1",  
           colNames=TRUE, 
           rowNames=TRUE, 
           append=FALSE, 
           showNA=TRUE)

cat(file=stderr(), "Summary results saved", "\n")

# ------------------------------------------------------------------------------
# plot data

cell_plots <- plot_cell_measurements(cell_measure_filter,
                                     plots_distance,
                                     measureChannelCell,
                                     measureChannelOrganelle,
                                     plot_background_subtract)

# TODO: This throws and error since it contains empty values
detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                               merged_summary,
                                               plots_distance,
                                               measureChannelCell,
                                               measureChannelOrganelle,
                                               cal_distance_nucleus,
                                               norm_distance_nucleus,
                                               plot_background_subtract)

do.call(grid.arrange, cell_plots)
do.call(grid.arrange, detection_plots)





# ------------------------------------------------------------------------------
# plot intensity profiles 

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
  
  write.xlsx(file = paste0( result_path,  "_intensityRatio.xlsx", sep = ""), 
             intensity_ratio_results, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
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

