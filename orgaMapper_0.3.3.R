library(ggplot2)
library(gridExtra)
library(tidyverse)
library("openxlsx")
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
#
#      VERSION: 0.3.3
#      CREATED: 2020-07-28
#     REVISION: 2020-10-12
#
# ==============================================================================
# user defined parameters

# path to folder where the directories for the measurements are
directory = "/home/schmiedc/Desktop/Test/test_nd2/output/"

result_name = "Analysis_test"

# filter for feret's diameter
feret_lower = 0
feret_upper = 600

# determine range for plots
norm_distance_nucleus = 0.7

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

# get number of columns from datasets
orga_column = ncol(organelle_distance);
cell_column = ncol(cell_measure);

# background subtraction
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
plot <- ggplot(norm_list2, aes(x = x, y = peak_norm, group = name, color = name)) + 
  geom_line() +
  xlab("Normalized distance from Nucleus") +
  ylab("Lysosome density (peak norm)") +
  ggtitle("Peak normalized density plots") +
  scale_x_continuous(expand = c(0, 0)) + # force start at 0
  scale_y_continuous(expand = c(0, 0)) + # force start at 0
  lineplot_theme()
  

print(plot)

ggsave(file=paste0(plots, .Platform$file.sep, "idensityPlot", ".pdf"), 
       width = 297, height = 210, units = "mm")

# ==============================================================================
# get value distances for each image
nameIntMeas = "intDistance.csv"
lysDist.filenames <- list.files(recursive=TRUE, pattern=nameIntMeas,full.names = TRUE)
head(lysDist.filenames)
value.list <- list()

head(merge.table_intensity)

intMeasure <- read.csv(lysDist.filenames[1], header = TRUE)
head(intMeasure)

for (name in lysDist.filenames[1]) {
  
  intMeasure <- read.csv(name, header = TRUE)
  intColumn = ncol(intMeasure)
  
  print(">>Start 1")
  # merge cell measurements with intensity profiles
  merge.table_intensity <- merge(cell_measure_filter, 
                                 intMeasure, 
                                 by = c("Name", "Series", "Cell"))
  
  # distance normalization
  merge.table_intensity$DistanceNorm <- merge.table_intensity$DistanceCal / merge.table_intensity$Ferets
  
  # background subtraction for detection intensity
  merge.table_intensity$orgaIntBackSub <- merge.table_intensity$orgaInt - merge.table_intensity$MeanBackgroundOrga
  
  if (cell_column == 10 && intColumn == 9) {
    
    merge.table_intensity$measureIntBackSub <- merge.table_intensity$measureInt - merge.table_intensity$MeanBackgroundMeasure
    
  }
  
  head(merge.table_intensity)
  
  # calculate mean intensity over bins with fixed width between 0 and 0.1
  orga <- tapply(merge.table_intensity$orgaInt, 
                 cut(merge.table_intensity$DistanceCal, 
                     seq(0, 0.1, by=0.005)), mean)
  
  head(orga)
  
  # transform output of tapply to data frame
  orga1 <- as.data.frame(as.table(orga))
  colnames(orga1) <- c("bin", "mean_orga")
  
  
  
  print(">>Start 2")
  # cleanup data frame
  orga1$bin <- str_remove_all(orga1$bin, "[]\\(]")
  orga1$bin <- str_replace(orga1$bin, ",", "-")
  
  # calculate mean intensity over bins with fixed width between 0 and 0.1
  orga_norm <- tapply(merge.table_intensity$orgaIntBackSub, 
                      cut(merge.table_intensity$DistanceNorm, 
                          seq(0, 0.1, by=0.005)), mean)
  
  # transform output of tapply to data frame
  orga_norm1 <- as.data.frame(as.table(orga_norm))
  colnames(orga_norm1) <- c("bin", "mean_orga_norm")
  
  # cleanup data frame
  orga_norm1$bin <- str_remove_all(orga_norm1$bin, "[]\\(]")
  orga_norm1$bin <- str_replace(orga_norm1$bin, ",", "-")
  
  # remerge with name and other data 
  value_result <- merge (orga1, orga_norm1, by = "bin")
  
  if (cell_column == 10 && intColumn == 9) {
    
    print(">>Start 3")
    # calculate mean intensity over bins with fixed width between 0 and 0.1
    measure <- tapply(merge.table_intensity$measureInt, 
                      cut(merge.table_intensity$DistanceCal, 
                          seq(0, 0.1, by=0.005)), mean)
    
    # transform output of tapply to data frame
    measure1 <- as.data.frame(as.table(measure))
    colnames(measure1) <- c("bin", "mean_measure")
    
    # cleanup data frame
    measure1$bin <- str_remove_all(measure1$bin, "[]\\(]")
    measure1$bin <- str_replace(measure1$bin, ",", "-")
    
    print(">>Start 4")
    # calculate mean intensity over bins with fixed width between 0 and 0.1
    measure_norm <- tapply(merge.table_intensity$measureIntBackSub, 
                           cut(merge.table_intensity$DistanceNorm, 
                               seq(0, 0.1, by=0.005)), mean)
    
    # transform output of tapply to data frame
    measure_norm1 <- as.data.frame(as.table(measure_norm))
    colnames(measure_norm1) <- c("bin", "mean_measure_norm")
    
    # cleanup data frame
    measure_norm1$bin <- str_remove_all(measure_norm1$bin, "[]\\(]")
    measure_norm1$bin <- str_replace(measure_norm1$bin, ",", "-")
    
    # remerge with name and other data 
    measure_result <- merge(measure1, measure_norm1, by = "bin")
    value_result <- merge(value_result, measure_result, by = "bin")
    
  }
  
  table.info <- merge.table_intensity[1,1:12]
  value_result2 <- merge(table.info, value_result, by=NULL)
  
  value.list[[name]] <- value_result2
  
}

# binds collection of normalized density plots and binds them into one dataframe
value.listCollected <- do.call("rbind", value.list)
head(value.listCollected)

# ------------------------------------------------------------------------------
# merged_summary
write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
           value.listCollected, 
           sheetName="Sheet1",  
           col.names=TRUE, 
           row.names=TRUE, 
           append=FALSE, 
           showNA=TRUE)

# ------------------------------------------------------------------------------
head(value.listCollected)

# calculate mean per cell
summary.intensity <- value.listCollected %>% 
  group_by(Name, bin) %>% 
  summarise(across(mean_intensity, mean, na.rm =TRUE ))

head(summary.intensity)

plot <- ggplot(data=summary.intensity, aes(x=bin, y=mean_intensity, group=Name)) +
  geom_line(aes(color=Name)) +
  geom_point(aes(color=Name)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Normalized distance from Nucleus") +
  ylab("Fluorescent intensity (A.U.)") +
  ggtitle("Intensity profile")

print(plot)

ggsave(file=paste0(plots, .Platform$file.sep, "intensityProfile_orgaCh", ".pdf"), 
       width = 297, height = 210, units = "mm")


















# ==============================================================================
# Statistical analysis
# split data into control and treatment
controlName = 'HeLa_scr'
treatName = 'HeLa_siArl8b'

HeLa_scr <- subset(merged_summary, Name == controlName)
HeLa_siArl8b <- subset(merged_summary, Name == treatName)

# DistanceNorm.mean
## Shapiro-Wilk Normality Test - normality of univariante sample
with(merged_summary, shapiro.test(DistanceNorm.mean[Name == controlName]))
with(merged_summary, shapiro.test(DistanceNorm.mean[Name == treatName]))

## if data is normally distributed then perform F-test
## f-test checks if both distributions have similar variance

## not normaly distributed thus Wilcoxon Test
test <- wilcox.test(HeLa_scr$DistanceNorm.mean, HeLa_siArl8b$DistanceNorm.mean)
test$p.value

# DistanceCal.mean
with(merged_summary, shapiro.test(DistanceCal.mean[Name == controlName]))
with(merged_summary, shapiro.test(DistanceCal.mean[Name == treatName]))

test <- wilcox.test(HeLa_scr$DistanceCal.mean, HeLa_siArl8b$DistanceCal.mean)
test$p.value

# CellArea
with(merged_summary, shapiro.test(CellArea[Name == controlName]))
with(merged_summary, shapiro.test(CellArea[Name == treatName]))

test <- wilcox.test(HeLa_scr$CellArea, HeLa_siArl8b$CellArea)
test$p.value

# Ferets
with(merged_summary, shapiro.test(Ferets[Name == controlName]))
with(merged_summary, shapiro.test(Ferets[Name == treatName]))

test <- wilcox.test(HeLa_scr$Ferets, HeLa_siArl8b$Ferets)
test$p.value

# NumDetections
with(merged_summary, shapiro.test(NumDetections[Name == controlName]))
with(merged_summary, shapiro.test(NumDetections[Name == treatName]))

test <- wilcox.test(HeLa_scr$NumDetections, HeLa_siArl8b$NumDetections)
test$p.value

