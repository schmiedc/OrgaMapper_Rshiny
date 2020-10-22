setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")
source("process_data.R")

# ==============================================================================
# Params

# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/Test/test_nd2/2020-10-14_output/"

result_name = "Analysis_test"

# filter for feret's diameter
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

orga_column = ncol(organelle_distance)
cell_column = ncol(cell_measure)

cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                 feret_lower, 
                                                 feret_upper,
                                                 cell_column,
                                                 orga_column)

merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                  organelle_distance,
                                                  cell_column,
                                                  orga_column)

merged_summary <- create_summary_table(merge_cell_organelle,
                                       cell_measure_filter)
