library(ggplot2)
library(tidyverse)
library("openxlsx")
library(data.table)
library(gridExtra)
# ==============================================================================
#
#  DESCRIPTION: Analysis script for OrgaMapper
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES:
#
#      VERSION: 0.4.1
#      CREATED: 2020-07-28
#     REVISION: 2020-10-12
#
# ==============================================================================
# user defined parameters

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
# where to save the data
setwd(directory)
input_dir = getwd()
out_dir =  getwd()

result_path <- file.path(out_dir, result_name, fsep = .Platform$file.sep)

# plot dir
plots = "plots"
dir.create("plots", showWarnings = FALSE)

# ==============================================================================
# data processing cells and detection data
# ==============================================================================
# read the data
name_distance = "organelleDistance.csv"
name_cell_measure = "cellMeasurements.csv"

organelle_distance <- read.csv(name_distance, header = TRUE)
cell_measure <- read.csv(name_cell_measure, header = TRUE)

# ------------------------------------------------------------------------------
# deal with different name and series options

# TODO needs to default to something that is sensible if invalid
if (single_series) {
  
  organelle_distance$Series <- str_extract(organelle_distance$Name, series_regex)
  cell_measure$Series <- str_extract(cell_measure$Name, series_regex)
  
  organelle_distance$Name <- str_remove(organelle_distance$Name, series_regex)
  cell_measure$Name <- str_remove(cell_measure$Name, series_regex)
  
  # removes trailing underscore or hyphen
  organelle_distance$Name <- str_remove(organelle_distance$Name, "(_|-| )($)")
  cell_measure$Name <- str_remove(cell_measure$Name, "(_|-| )($)")
  
}

# ------------------------------------------------------------------------------
# get number of columns from datasets
orga_column = ncol(organelle_distance);
cell_column = ncol(cell_measure);

# ------------------------------------------------------------------------------
# data processing
cell_measure_filter <- subset(cell_measure, 
                             Ferets >= feret_lower & Ferets <= feret_upper)

# background subtraction for mean organelle intensity per cell
cell_measure_filter$MeanOrgaBackSub <- 
  cell_measure_filter$MeanValueOrga - cell_measure_filter$MeanBackgroundOrga

if (cell_column == 10 && orga_column == 10) {
  
  cell_measure_filter$MeanMeasureBackSub <- 
    cell_measure_filter$MeanValueMeasure - cell_measure_filter$MeanBackgroundMeasure
  
}

merge_cell_organelle <- merge(cell_measure_filter, 
                              organelle_distance, 
                              by = c("Name", "Series", "Cell"))

# background subtraction for detection intensity
merge_cell_organelle$PeakDetectBackSub <- 
  merge_cell_organelle$PeakDetectionInt - merge_cell_organelle$MeanBackgroundOrga

if (cell_column == 10 && orga_column == 10) {
  
  merge_cell_organelle$PeakMeasureBackSub <- 
    merge_cell_organelle$PeakMeasureInt - merge_cell_organelle$MeanBackgroundMeasure
  
}

merge_cell_organelle$DistanceNorm <- 
  merge_cell_organelle$DistanceCal / merge_cell_organelle$Ferets


summary_table <- merge_cell_organelle %>% 
  group_by(Name, Series, Cell) %>% 
  summarise(across(DistanceRaw:DistanceNorm, mean, na.rm =TRUE )) %>% 
  rename_at(vars(-Name, -Series, -Cell),function(x) paste0(x,".mean"))

merged_summary <- merge(cell_measure_filter, 
                        summary_table, by = c("Name", "Series", "Cell"))

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
# ==============================================================================
# plot detection and cell data
# ==============================================================================
boxplot_theme <- function() {
  theme_bw()
  
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # x axis ticks
        axis.text.x = element_text(color = "grey20", 
                                   size = 10, 
                                   angle = 45, 
                                   hjust = .5, 
                                   vjust = .5, 
                                   face = "plain"),
        # y axis ticks
        axis.text.y = element_text(color = "grey20", 
                                   size = 10, 
                                   angle = 0, 
                                   hjust = 1, 
                                   vjust = 0, 
                                   face = "plain"),
        # x axis labels
        axis.title.x = element_text(color = "grey20", 
                                    size = 12, 
                                    angle = 0, 
                                    hjust = .5, 
                                    vjust = 0, 
                                    face = "plain"),
        # y axis labels
        axis.title.y = element_text(color = "grey20", 
                                    size = 12, 
                                    angle = 90, 
                                    hjust = .5, 
                                    vjust = .5, 
                                    face = "plain"),
        # title
        title = element_text(color = "grey20", 
                             size = 14, 
                             angle = 0, 
                             hjust = 0, 
                             vjust = 1, 
                             face = "plain")
  )
}

lineplot_theme <- function() {
  theme_bw()
  
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # x axis ticks
    axis.text.x = element_text(color = "grey20", 
                               size = 10, 
                               angle = 0, 
                               hjust = .5, 
                               vjust = .5, 
                               face = "plain"),
    # y axis ticks
    axis.text.y = element_text(color = "grey20", 
                               size = 10, 
                               angle = 0, 
                               hjust = 1, 
                               vjust = 0, 
                               face = "plain"),
    # x axis labels
    axis.title.x = element_text(color = "grey20", 
                                size = 12, 
                                angle = 0, 
                                hjust = .5, 
                                vjust = 0, 
                                face = "plain"),
    # y axis labels
    axis.title.y = element_text(color = "grey20", 
                                size = 12, 
                                angle = 90, 
                                hjust = .5, 
                                vjust = .5, 
                                face = "plain"),
    # title
    title = element_text(color = "grey20", 
                         size = 14, 
                         angle = 0, 
                         hjust = 0, 
                         vjust = 1, 
                         face = "plain")
  )
}
# ------------------------------------------------------------------------------
# plots area, number of detections and mean value
organelle_intensity_cell = ""

if (plot_background_subtract) {
  
  organelle_intensity_cell = "MeanOrgaBackSub"
  
} else {
  
  organelle_intensity_cell = "MeanValueOrga"
  
}

measure1 <- list("Ferets", 
                 "CellArea", 
                 "NumDetections",
                 organelle_intensity_cell)

measure1_title <- list("Average Feret's diameter", 
                       "Average cell area", 
                       "Number of detections per cell", 
                       "Average intensity in cell (detection channel)")

measure1_label <- list("ferets diameter (µm)", 
                       "cell area (µm²)", 
                       "average count", 
                       "fluorescent intensity (A.U.)")

measure1_file <- list("feretPlot", 
                      "cellArea", 
                      "numDetections", 
                      "intPerCellDetection")

cell_measure_long <- cell_measure_filter %>% 
  pivot_longer(cols=Ferets:MeanOrgaBackSub,values_to = "measurement" )

plot_list_cell <- list()

for (index in seq_along(measure1)) {
  
  dataSubset <- subset(cell_measure_long, cell_measure_long$name==measure1[index])
  
  plot_cell <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    ggtitle(measure1_title[index]) + 
    xlab("Treatment") +
    ylab( measure1_label[index] ) + 
    boxplot_theme()

  plot_list_cell[[index]] <- plot_cell 
  
  ggsave(plot = plot_cell,
         file=paste0(plots, .Platform$file.sep, measure1_file[index], ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
}

# ------------------------------------------------------------------------------
# lysosome density plots

plot_list_detection <- list()

name_count <- as.data.frame(table(merge_cell_organelle$Name))
detect_list <- list()

# goes through each experiment and calculates lysosome density
# then peak normalizes the lysosome density
# collects these normalized density plots in detect_list
for (name in name_count$Var1){
  
  data_per_name <- subset(merge_cell_organelle, Name == name)
  density_per_name <- density(data_per_name$DistanceNorm, 
                              bw = "nrd0", 
                              n = 512, 
                              from = 0, 
                              to = norm_distance_nucleus)
  
  data_frame <- data.frame(density_per_name$x)
  data_frame$y <- density_per_name$y
  colnames(data_frame)[1] <-  "x"
  
  # peak normalisation
  max = max(data_frame$y, na.rm = FALSE)
  data_frame$peak_norm <- sapply(data_frame$y, function(x){x /  max})
  detect_list[[name]] <- data_frame
  
}

# binds collection of normalized density plots and binds them into one dataframe
norm_list <- do.call("rbind", detect_list)
norm_list1 <- tibble::rownames_to_column(norm_list, "nameindex")
norm_list1_indices <- str_split_fixed(norm_list1$nameindex, "\\.", 2)
norm_list2 <- cbind(norm_list1_indices, norm_list1)
colnames(norm_list2)[1] <- "name"
colnames(norm_list2)[2] <- "index"

# Plot Lysosome density vs normalized distance from Nucleus
# density plots without peak normalized data
plot_density_raw <- ggplot(norm_list2, aes(x = x, 
                                           y = y, 
                                           group = name, 
                                           color = name)) + 
  geom_line() +
  xlab("Normalized distance from Nucleus") +
  ylab("Lysosome density (peak norm)") +
  ggtitle("Raw density plots") +
  scale_x_continuous(expand = c(0, 0)) + # force start at 0
  scale_y_continuous(expand = c(0, 0)) + # force start at 0
  lineplot_theme()

ggsave(plot = plot_density_raw,
       file=paste0(plots, .Platform$file.sep, "raw_densityPlot", ".pdf"), 
       width = 297, 
       height = 210, 
       units = "mm")

plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_density_raw

# Plot Lysosome density vs normalized distance from Nucleus
# density plots with peak normalized data
plot_density <- ggplot(norm_list2, aes(x = x, 
                                       y = peak_norm, 
                                       group = name, 
                                       color = name)) + 
  geom_line() +
  xlab("Normalized distance from Nucleus") +
  ylab("Lysosome density (peak norm)") +
  ggtitle("Peak normalized density plots") +
  scale_x_continuous(expand = c(0, 0)) + # force start at 0
  scale_y_continuous(expand = c(0, 0)) + # force start at 0
  lineplot_theme()

ggsave(plot = plot_density,
       file=paste0(plots, .Platform$file.sep, "densityPlot", ".pdf"), 
       width = 297, 
       height = 210, 
       units = "mm")

plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_density

# ------------------------------------------------------------------------------
# plot distance and detection intensity
head(summary_table)
organelle_intensity_cell = ""

if (plot_background_subtract) {
  
  organelle_intensity_peak = "PeakDetectBackSub.mean"
  
} else {
  
  organelle_intensity_peak = "PeakDetectionInt.mean"
  
}

measure2 <- list("DistanceCal.mean", 
                 "DistanceNorm.mean",
                 organelle_intensity_peak)

measure2_title <- list("Average distance from nucleus", 
                       "Average normalized distance from nucleus",
                       "Average peak detection intensity  (detection channel)")

measure2_label <- list("distance (µm)", 
                       "normalized distance",
                       "fluorescent intensity (A.U.)")

measure2_file <- list("distance", 
                      "distanceNorm",
                      "intPerDetectionDetectionChannel")

summary_long <- summary_table %>% 
  pivot_longer(cols=DistanceRaw.mean:DistanceNorm.mean, 
               values_to = "measurement" )


for (index in seq_along(measure2)) {
  
  dataSubset <- subset(summary_long, summary_long$name==measure2[index])
  
  plot_detection <- distancePlot <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    ggtitle(measure2_title[index]) + 
    xlab("Treatment") +
    ylab(measure2_label[index]) +
    boxplot_theme()
  
  plot_list_detection[[length(plot_list_detection) + 1]] <- plot_detection
  
  ggsave(plot = plot_detection, 
         file=paste0(plots, .Platform$file.sep, measure2_file[index], ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
}


# ------------------------------------------------------------------------------
# plot measure channel

measure_intensity_cell = ""
measure_intensity_peak = ""

if (cell_column == 10 && orga_column == 10) {
  
  if (plot_background_subtract) {
    
    measure_intensity_cell = "MeanMeasureBackSub"
    measure_intensity_peak = "PeakMeasureBackSub.mean"
    
  } else {
    
    measure_intensity_cell = "MeanValueMeasure"
    measure_intensity_peak = "PeakMeasureInt.mean"
    
  }
  
  head(cell_measure_filter)
  
  cell_measure_filter_new <- cell_measure_filter[c("Name", measure_intensity_cell)]
  colnames(cell_measure_filter_new)[2] <- "measure"
  
  plot_cell_measure <- ggplot(cell_measure_filter_new, aes(x=Name, y=measure)) +
    ggtitle("Average intensity in cell (measure channel)") + 
    xlab("Treatment") +
    ylab("fluorescent intensity (A.U.)") +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    boxplot_theme()
  
  plot_list_cell[[length(plot_list_cell)  + 1]] <- plot_cell_measure
  
  ggsave(plot = plot_cell_measure,
         file=paste0(plots, .Platform$file.sep, "intPerCellMeasure", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  summary_table_new <- summary_table[c("Name", measure_intensity_peak)]
  colnames(summary_table_new)[2] <- "measure"
  
  plot_peak_measure <- ggplot(summary_table_new, aes(x=Name, y=measure)) +
    ggtitle("Average peak detection intensity (measure channel)") + 
    xlab("Treatment") +
    ylab("fluorescent intensity (A.U.)") +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    boxplot_theme()
  
  plot_list_detection[[length(plot_list_detection)  + 1]] <- plot_peak_measure
  
  ggsave(plot = plot_peak_measure,
         file=paste0(plots, 
                     .Platform$file.sep, 
                     "intPerDetectionMeasureChannel", 
                     ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
}

# ==============================================================================
# print plots to grid
# ==============================================================================
do.call(grid.arrange, plot_list_cell)
do.call(grid.arrange, plot_list_detection)


# ==============================================================================
# intensity profiles
# ==============================================================================
if (analyze_signal_profiles) {
  
  # get value distances for each image
  name_value_measure = "intDistance.csv"
  value_filenames <- list.files(recursive=TRUE, 
                                  pattern=name_value_measure,full.names = TRUE)
  
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
  
  value_list_norm <- list()
  value_list <- list()
  
  for (name in value_filenames) {
    
    value_measure  <- read.csv(name, header = TRUE)
    
    # TODO needs to default to something that is sensible if invalid
    if (single_series) {
      
      value_measure$Series <- str_extract(value_measure$Name, series_regex)
      
      value_measure$Name <- str_remove(value_measure$Name, series_regex)
      
      # removes trailing underscore or hypen
      value_measure$Name <- str_remove(value_measure$Name, "(_|-| )($)")
      
    }
    
    value_column = ncol(value_measure )
    
    # merge cell measurements with intensity profiles
    merge_table_value <- merge(cell_measure_filter, 
                                   value_measure , 
                                   by = c("Name", "Series", "Cell"))
    
    # distance normalization
    merge_table_value$DistanceNorm <- merge_table_value$DistanceCal / merge_table_value$Ferets
    
    # background subtraction for detection intensity
    merge_table_value$orgaIntBackSub <- merge_table_value$orgaInt - merge_table_value$MeanBackgroundOrga
    
    if (cell_column == 10 && value_column == 9) {
      
      merge_table_value$measureIntBackSub <- merge_table_value$measureInt - merge_table_value$MeanBackgroundMeasure
      
    }
    
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
    
    if (cell_column == 10 && value_column == 9) {
      
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
  
  # ------------------------------------------------------------------------------
  # save intensity profiles
  # ------------------------------------------------------------------------------
  write.xlsx(file = paste0( result_path,  "_intensityProfile_norm.xlsx", sep = ""), 
             value_list_norm_coll, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
             value_list_coll, 
             sheetName="Sheet1",  
             col.names=TRUE, 
             row.names=TRUE, 
             append=FALSE, 
             showNA=TRUE)
  
  # ----------------------------------------------------------------------------
  # plot intensity profiles organelle channel
  # ----------------------------------------------------------------------------
  intensity_profile_plotlist <- list()
  
  summary_value_norm <- data.frame(Date=as.Date(character()),
                                   File=character(), 
                                   User=character(), 
                                   stringsAsFactors=FALSE) 
  
  summary_value <- data.frame(Date=as.Date(character()),
                               File=character(), 
                               User=character(), 
                               stringsAsFactors=FALSE) 
  
  if (plot_background_subtract) {
    
    summary_value_norm <- value_list_norm_coll %>% 
      group_by(Name, bin) %>% 
      summarise(across(backsub_bin_orga_norm, mean, na.rm =TRUE ))
    
    summary_value <- value_list_coll %>% 
      group_by(Name, bin, row) %>% 
      summarise(across(backsub_bin_orga, mean, na.rm =TRUE ))
    
  } else {
    
    summary_value_norm <- value_list_norm_coll %>% 
      group_by(Name, bin) %>% 
      summarise(across(raw_bin_orga_norm, mean, na.rm =TRUE ))
    
    summary_value <- value_list_coll %>% 
      group_by(Name, bin, row) %>% 
      summarise(across(raw_bin_orga, mean, na.rm =TRUE ))
    
  }
  
  names(summary_value_norm)[3] <- "value"
  names(summary_value)[4] <- "value"
  
  # calculate mean per cell for raw data
  plot_profile <- ggplot(data=summary_value, aes(x=reorder(bin,row), y=value, group=Name)) +
    geom_line(aes(color=Name)) +
    geom_point(aes(color=Name)) +
    lineplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    ggtitle("Intensity profile organelle Channel") +
    xlab("Distance from Nucleus (µm)") 
  
  ggsave(plot = plot_profile,
         file=paste0(plots, .Platform$file.sep, "intensityProfile_orgaCh", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  intensity_profile_plotlist[[length(intensity_profile_plotlist) + 1 ]] <- plot_profile
  # calculate mean per cell for normalized data

  plot_profile_norm <- ggplot(data=summary_value_norm, aes(x=bin, y=value, group=Name)) +
    geom_line(aes(color=Name)) +
    geom_point(aes(color=Name)) +
    lineplot_theme() + 
    aes(x = fct_inorder(bin)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Fluorescent intensity (A.U.)") +
    xlab("Normalized distance from Nucleus") +
    ggtitle("Intensity profile organelle Channel")

  ggsave(plot = plot_profile_norm,
         file=paste0(plots, .Platform$file.sep, "NormIntensityProfile_orgaCh", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  intensity_profile_plotlist[[length(intensity_profile_plotlist) + 1 ]] <- plot_profile_norm
  
  # peak normalisation
  name_count_value <- as.data.frame(table(summary_value_norm$Name))
  value_list <- list()
  
  for (name in name_count_value$Var1){
  
    data_per_name <- subset(summary_value_norm, Name == name)
    max_value_profiles = max(data_per_name$value, na.rm = TRUE)
    data_per_name$peak_norm <- sapply(data_per_name$value, function(x){x /  max_value_profiles})
    
    value_list[[name]] <- data_per_name
  }
  # binds collection of normalized value plots and binds them into one dataframe
  norm_list_value <- do.call("rbind", value_list)
  
  plot_profile_peak <- ggplot(data=norm_list_value, aes(x=bin, y=peak_norm, group=Name)) +
    geom_line(aes(color=Name)) +
    geom_point(aes(color=Name)) +
    lineplot_theme() + 
    aes(x = fct_inorder(bin)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Normalized fluorescent intensity") +
    xlab("Normalized distance from Nucleus") +
    ggtitle("Peak normalized intensity profile organelle Channel")

  ggsave(plot = plot_profile_peak,
         file=paste0(plots, .Platform$file.sep, "PeakNormIntensityProfile_orgaCh", ".pdf"), 
         width = 297, 
         height = 210, 
         units = "mm")
  
  intensity_profile_plotlist[[length(intensity_profile_plotlist) + 1 ]] <- plot_profile_peak
  
  do.call(grid.arrange, intensity_profile_plotlist)
  
  head(value_list_norm_coll)
  head(value_list_coll)
  # ----------------------------------------------------------------------------
  # plot intensity profiles measure channel
  # ----------------------------------------------------------------------------
  if (cell_column == 10 && value_column == 9) {
    
    intensity_profile_plotlist_meas <- list()
    
    summary_value_norm_meas <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 
    
    summary_value_meas <- data.frame(Date=as.Date(character()),
                                     File=character(), 
                                     User=character(), 
                                     stringsAsFactors=FALSE) 
    
    if (plot_background_subtract) {
      
      summary_value_norm_meas <- value_list_norm_coll %>% 
        group_by(Name, bin) %>% 
        summarise(across(backsub_bin_measure_norm, mean, na.rm =TRUE ))
      
      summary_value_meas <- value_list_coll %>% 
        group_by(Name, bin, row) %>% 
        summarise(across(backsub_bin_measure, mean, na.rm =TRUE ))
      
    } else {
      
      summary_value_norm_meas <- value_list_norm_coll %>% 
        group_by(Name, bin) %>% 
        summarise(across(raw_bin_measure_norm, mean, na.rm =TRUE ))
      
      summary_value_meas <- value_list_coll %>% 
        group_by(Name, bin, row) %>% 
        summarise(across(raw_bin_measure, mean, na.rm =TRUE ))
      
    }
    
    names(summary_value_norm_meas)[3] <- "value"
    names(summary_value_meas)[4] <- "value"
    
    # calculate mean per cell for raw data
    plot_profile_meas <- ggplot(data=summary_value_meas, aes(x=reorder(bin,row), y=value, group=Name)) +
      geom_line(aes(color=Name)) +
      geom_point(aes(color=Name)) +
      lineplot_theme() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("Fluorescent intensity (A.U.)") +
      ggtitle("Intensity profile measure Channel") +
      xlab("Distance from Nucleus (µm)") 
    
    ggsave(plot = plot_profile_meas,
           file=paste0(plots, .Platform$file.sep, "intensityProfile_measCh", ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
    intensity_profile_plotlist_meas[[length(intensity_profile_plotlist_meas) + 1 ]] <- plot_profile_meas
    # calculate mean per cell for normalized data
    
    plot_profile_norm_meas <- ggplot(data=summary_value_norm_meas, aes(x=bin, y=value, group=Name)) +
      geom_line(aes(color=Name)) +
      geom_point(aes(color=Name)) +
      lineplot_theme() + 
      aes(x = fct_inorder(bin)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("Fluorescent intensity (A.U.)") +
      xlab("Normalized distance from Nucleus") +
      ggtitle("Intensity profile measure Channel")
    
    ggsave(plot = plot_profile_norm_meas,
           file=paste0(plots, .Platform$file.sep, "NormIntensityProfile_measCh", ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
    intensity_profile_plotlist_meas[[length(intensity_profile_plotlist_meas) + 1 ]] <- plot_profile_norm_meas
    
    # peak normalisation
    name_count_value <- as.data.frame(table(summary_value_norm_meas$Name))
    value_list <- list()
    
    for (name in name_count_value$Var1){
      
      data_per_name <- subset(summary_value_norm_meas, Name == name)
      max_value_profiles = max(data_per_name$value, na.rm = TRUE)
      data_per_name$peak_norm <- sapply(data_per_name$value, function(x){x /  max_value_profiles})
      
      value_list[[name]] <- data_per_name
    }
    # binds collection of normalized value plots and binds them into one dataframe
    norm_list_value <- do.call("rbind", value_list)
    
    plot_profile_peak_meas <- ggplot(data=norm_list_value, aes(x=bin, y=peak_norm, group=Name)) +
      geom_line(aes(color=Name)) +
      geom_point(aes(color=Name)) +
      lineplot_theme() + 
      aes(x = fct_inorder(bin)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("Normalized fluorescent intensity") +
      xlab("Normalized distance from Nucleus") +
      ggtitle("Peak normalized intensity profile measure Channel")
    
    ggsave(plot = plot_profile_peak_meas,
           file=paste0(plots, .Platform$file.sep, "PeakNormIntensityProfile_measCh", ".pdf"), 
           width = 297, 
           height = 210, 
           units = "mm")
    
    intensity_profile_plotlist_meas[[length(intensity_profile_plotlist_meas) + 1 ]] <- plot_profile_peak_meas
    
    do.call(grid.arrange, intensity_profile_plotlist_meas)
    
  }

}
