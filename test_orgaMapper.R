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
directory = "/data2/shared_data/OrgaMapper_Data/size_MTM1KOvsWT/output/"

result_name = "Analysis_test"

# filter for feret's diameter
feret_filter = TRUE
feret_lower = 0
feret_upper = 600

# determine range for plots
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
plots_dir <- file.path(out_dir, "plots", fsep = .Platform$file.sep)
dir.create(plots_dir, showWarnings = FALSE)

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
# save processed data
write.xlsx(file = paste0( result_path, "_detection.xlsx", sep = ""), 
           merge_cell_organelle, 
           sheetName="Sheet1",  
           col.names=TRUE, 
           row.names=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# merged_summary
write.xlsx(file = paste0( result_path,  "_cell.xlsx", sep = ""), 
           merged_summary, 
           sheetName="Sheet1",  
           col.names=TRUE, 
           row.names=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# ------------------------------------------------------------------------------
# plot data
cell_plots <- plot_cell_measurements(cell_measure_filter,
                                     plots_dir,
                                     cell_column,
                                     orga_column,
                                     plot_background_subtract)

detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                               merged_summary,
                                               plots_dir,
                                               cell_column,
                                               orga_column,
                                               norm_distance_nucleus,
                                               plot_background_subtract)

do.call(grid.arrange, cell_plots)
do.call(grid.arrange, detection_plots)

# if (analyze_signal_profiles) {
  
  name_value_measure = "intDistance.csv"
  
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

  # ------------------------------------------------------------------------------
  write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
             value_list, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  write.xlsx(file = paste0( result_path,  "_intensityProfile_norm.xlsx", sep = ""), 
             value_list_norm, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  profile_plot <- plot_profiles(value_list, 
                                value_list_norm, 
                                "organelle", 
                                plots_dir, 
                                plot_background_subtract)
  
  if (cell_column == 10 && orga_column == 10) {
    
    measure_profiles <- plot_profiles(value_list, 
                                      value_list_norm, 
                                      "measure", 
                                      plots_dir, 
                                      plot_background_subtract)
    
    do.call(grid.arrange, measure_profiles)
    
  }
  
  do.call(grid.arrange, profile_plot)
  
# }

