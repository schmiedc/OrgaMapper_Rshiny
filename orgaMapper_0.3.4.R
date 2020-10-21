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
#      VERSION: 0.3.4
#      CREATED: 2020-07-28
#     REVISION: 2020-10-12
#
# ==============================================================================
# user defined parameters

# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/OrgaMapper_useCases/SingelSeries_tiff/output/"

result_name = "Analysis_test"

# filter for feret's diameter
feret_lower = 0
feret_upper = 600

# determine range for plots
norm_distance_nucleus = 0.7

# TODO if file contains series number or the already present column
# needs to default to something sensible if not possible
single_series = TRUE
series_regex = "(?<=_)\\d*($)"

# TODO apply background subtraction for plots
background_subtract_plots = FALSE

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
# read the data
name_distance = "organelleDistance.csv"
name_cell_measure = "cellMeasurements.csv"

organelle_distance <- read.csv(name_distance, header = TRUE)
cell_measure <- read.csv(name_cell_measure, header = TRUE)

# ==============================================================================
# deal with different name and series options

# TODO needs to default to something that is sensible if invalid
if (single_series) {
  
  organelle_distance$Series <- str_extract(organelle_distance$Name, series_regex)
  cell_measure$Series <- str_extract(cell_measure$Name, series_regex)
  
  organelle_distance$Name <- str_remove(organelle_distance$Name, series_regex)
  cell_measure$Name <- str_remove(cell_measure$Name, series_regex)
  
  # removes trailing underscore or hypen
  organelle_distance$Name <- str_remove(organelle_distance$Name, "(_|-| )($)")
  cell_measure$Name <- str_remove(cell_measure$Name, "(_|-| )($)")
  
}

# ==============================================================================
# get number of columns from datasets
orga_column = ncol(organelle_distance);
cell_column = ncol(cell_measure);

# ==============================================================================
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

# ==============================================================================
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
                                   size = 16, 
                                   angle = 45, 
                                   hjust = .5, 
                                   vjust = .5, 
                                   face = "plain"),
        # y axis ticks
        axis.text.y = element_text(color = "grey20", 
                                   size = 16, 
                                   angle = 0, 
                                   hjust = 1, 
                                   vjust = 0, 
                                   face = "plain"),
        # x axis labels
        axis.title.x = element_text(color = "grey20", 
                                    size = 22, 
                                    angle = 0, 
                                    hjust = .5, 
                                    vjust = 0, 
                                    face = "plain"),
        # y axis labels
        axis.title.y = element_text(color = "grey20", 
                                    size = 22, 
                                    angle = 90, 
                                    hjust = .5, 
                                    vjust = .5, 
                                    face = "plain"),
        # title
        title = element_text(color = "grey20", 
                             size = 25, 
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
                               size = 16, 
                               angle = 0, 
                               hjust = .5, 
                               vjust = .5, 
                               face = "plain"),
    # y axis ticks
    axis.text.y = element_text(color = "grey20", 
                               size = 16, 
                               angle = 0, 
                               hjust = 1, 
                               vjust = 0, 
                               face = "plain"),
    # x axis labels
    axis.title.x = element_text(color = "grey20", 
                                size = 22, 
                                angle = 0, 
                                hjust = .5, 
                                vjust = 0, 
                                face = "plain"),
    # y axis labels
    axis.title.y = element_text(color = "grey20", 
                                size = 22, 
                                angle = 90, 
                                hjust = .5, 
                                vjust = .5, 
                                face = "plain"),
    # title
    title = element_text(color = "grey20", 
                         size = 25, 
                         angle = 0, 
                         hjust = 0, 
                         vjust = 1, 
                         face = "plain")
  )
}
# ==============================================================================
# plots area, number of detections and mean value
measure1 <- list("Ferets", 
                 "CellArea", 
                 "NumDetections", 
                 "MeanOrgaBackSub")

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

for (index in seq_along(measure1)) {
  
  dataSubset <- subset(cell_measure_long, cell_measure_long$name==measure1[index])
  
  plot <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    ggtitle(measure1_title[index]) + 
    xlab("Treatment") +
    ylab( measure1_label[index] ) + 
    boxplot_theme()
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, measure1_file[index], ".pdf"), 
         width = 297, height = 210, units = "mm")
  
}

# ------------------------------------------------------------------------------
# plot distance and detection intensity
measure2 <- list("PeakDetectBackSub.mean", 
                 "DistanceCal.mean", 
                 "DistanceNorm.mean")

measure2_title <- list("Average peak detection intensity  (detection channel)", 
                       "Average distance from nucleus", 
                       "Average normalized distance from nucleus")

measure2_label <- list("fluorescent intensity (A.U.)", 
                       "distance (µm)", 
                       "normalized distance")

measure2_file <- list("intPerDetectionDetectionChannel", 
                      "distance", 
                      "distanceNorm")

summary_long <- summary_table %>% 
  pivot_longer(cols=DistanceRaw.mean:DistanceNorm.mean, 
               values_to = "measurement" )

for (index in seq_along(measure2)) {
  
  dataSubset <- subset(summary_long, summary_long$name==measure2[index])
  
  plot <- distancePlot <- ggplot(dataSubset, aes(x=Name, y=measurement)) +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    ggtitle(measure2_title[index]) + 
    xlab("Treatment") +
    ylab(measure2_label[index]) +
    boxplot_theme()
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, measure2_file[index], ".pdf"), 
         width = 297, height = 210, units = "mm")
  
}

# ------------------------------------------------------------------------------
# plot measure channel
if (cell_column == 10 && orga_column == 10) {
  
  plot <- ggplot(cell_measure_filter, aes(x=Name, y=MeanMeasureBackSub)) +
    ggtitle("Average intensity in cell (measure channel)") + 
    xlab("Treatment") +
    ylab("fluorescent intensity (A.U.)") +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    boxplot_theme()
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, "intPerCellMeasure", ".pdf"), 
         width = 297, height = 210, units = "mm")
  
  plot <- ggplot(summary_table, aes(x=Name, y=PeakMeasureBackSub.mean)) +
    ggtitle("Average peak detection intensity (measure channel)") + 
    xlab("Treatment") +
    ylab("fluorescent intensity (A.U.)") +
    geom_boxplot(outlier.size = 0, outlier.shape = 1) + 
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(width = 0.1) +
    boxplot_theme()
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, "intPerDetectionMeasureChannel", ".pdf"), 
         width = 297, height = 210, units = "mm")
  
}

# ------------------------------------------------------------------------------
# lysosome density plots
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
  colnames(data_frame)[1] <- "x"
  
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
# density plots with peak normalized data
plot <- ggplot(norm_list2, aes(x = x, 
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
  

print(plot)

ggsave(file=paste0(plots, .Platform$file.sep, "densityPlot", ".pdf"), 
       width = 297, height = 210, units = "mm")

# ==============================================================================
# intensity profiles
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
                                                  "mean_orga_norm",
                                                  bin_width_norm,
                                                  upper_limit_norm)
    
    value_result <- bin_distance_values(merge_table_value$DistanceCal, 
                                             merge_table_value$orgaIntBackSub, 
                                             "mean_orga",
                                             bin_width,
                                             upper_limit)
    
    if (cell_column == 10 && value_column == 9) {
      
      binned_measure_value_norm <- bin_distance_values(merge_table_value$DistanceNorm, 
                                                    merge_table_value$measureIntBackSub, 
                                                    "mean_measure_norm",
                                                    bin_width_norm,
                                                    upper_limit_norm)
      
      binned_measure_value <- bin_distance_values(merge_table_value$DistanceCal, 
                                               merge_table_value$measureIntBackSub, 
                                               "mean_measure",
                                               bin_width,
                                               upper_limit)
      
      
      value_result_norm <- merge(value_result_norm, binned_measure_value_norm, by = c("bin","row"))
      value_result <- merge(value_result, binned_measure_value, by = c("bin","row"))
      
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
  
  # ------------------------------------------------------------------------------
  # calculate mean per cell for normalized data
  summary.intensity <- value_list_norm_coll %>% 
    group_by(Name, bin) %>% 
    summarise(across(mean_orga_norm, mean, na.rm =TRUE ))
  
  plot <- ggplot(data=summary.intensity, aes(x=bin, y=mean_orga_norm, group=Name)) +
    geom_line(aes(color=Name)) +
    geom_point(aes(color=Name)) +
    xlab("Normalized distance from Nucleus") +
    ylab("Fluorescent intensity (A.U.)") +
    ggtitle("Intensity profile organelle Channel") +
    lineplot_theme() + 
    aes(x = fct_inorder(bin)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, "NormIntensityProfile_orgaCh", ".pdf"), 
         width = 297, height = 210, units = "mm")
  
  
  # calculate mean per cell for normalized data
  summary.intensity <- value_list_coll %>% 
    group_by(Name, bin, row) %>% 
    summarise(across(mean_orga, mean, na.rm =TRUE ))
  
  plot <- ggplot(data=summary.intensity, aes(x=reorder(bin,row), y=mean_orga, group=Name)) +
    geom_line(aes(color=Name)) +
    geom_point(aes(color=Name)) +
    ylab("Fluorescent intensity (A.U.)") +
    ggtitle("Intensity profile organelle Channel") +
    lineplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("Distance from Nucleus (µm)") 
  
  print(plot)
  
  ggsave(file=paste0(plots, .Platform$file.sep, "intensityProfile_orgaCh", ".pdf"), 
         width = 297, height = 210, units = "mm")
  
}
