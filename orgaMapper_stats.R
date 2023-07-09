library("openxlsx")

setwd("/home/schmiedc/FMP_Docs/Projects/OrgaMapper/Publication_orgaMapper/OrgaMapper_Data/size_MTM1KOvsWT/output_test/")

file <- read.xlsx("Analysis_test_cell.xlsx")
head(file)

# ==============================================================================
# Statistical analysis
# split data into control and treatment
unique(file$identifier)

controlName = 'wt_HBSS'
treatName = 'ko_HBSS'

control <- subset(file, identifier == controlName)
treatment <- subset(file, identifier == treatName)

# ==============================================================================
# General parameters for cells
# CellArea
with(file, shapiro.test(cell_area[identifier == controlName]))
with(file, shapiro.test(cell_area[identifier == treatName]))

test <- wilcox.test(control$cell_area, treatment$cell_area)
test$p.value

# Ferets
with(file, shapiro.test(ferets[identifier == controlName]))
with(file, shapiro.test(ferets[identifier == treatName]))

test <- wilcox.test(control$ferets, treatment$ferets)
test$p.value

# orga_numberOfDetections
with(file, shapiro.test(orga_numberOfDetections[identifier == controlName]))
with(file, shapiro.test(orga_numberOfDetections[identifier == treatName]))

test <- wilcox.test(control$orga_numberOfDetections, 
                    treatment$orga_numberOfDetections)
test$p.value

# average intensity
with(file, shapiro.test(orga_intensity_backsub[identifier == controlName]))
with(file, shapiro.test(orga_intensity_backsub[identifier == treatName]))

test <- wilcox.test(control$orga_intensity_backsub, 
                    treatment$orga_intensity_backsub)
test$p.value

# ==============================================================================
# orga_meanDistance_calibrated
## Shapiro-Wilk Normality Test - normality of univariante sample
with(file, shapiro.test(orga_meanDistance_calibrated[identifier == controlName]))
with(file, shapiro.test(orga_meanDistance_calibrated[identifier == treatName]))

## if data is normally distributed then perform F-test
## f-test checks if both distributions have similar variance

## not normaly distributed thus Wilcoxon Test
test <- wilcox.test(control$orga_meanDistance_calibrated, 
                    treatment$orga_meanDistance_calibrated)
test$p.value

# orga_meanDistance_normalized
with(file, shapiro.test(orga_meanDistance_normalized[identifier == controlName]))
with(file, shapiro.test(orga_meanDistance_normalized[identifier == treatName]))

test <- wilcox.test(control$orga_meanDistance_normalized, 
                    treatment$orga_meanDistance_normalized)
test$p.value

# orga_intensityOnDetection_backsub
with(file, shapiro.test(orga_intensityOnDetection_backsub[identifier == controlName]))
with(file, shapiro.test(orga_intensityOnDetection_backsub[identifier == treatName]))

test <- wilcox.test(control$orga_intensityOnDetection_backsub, 
                    treatment$orga_intensityOnDetection_backsub)
test$p.value

# ==============================================================================
ratio_file <- read.xlsx("Analysis_test_intensityRatio.xlsx", rowNames = TRUE)
head(ratio_file)

HeLa_scr_ratio <- subset(ratio_file, identifier == controlName)
HeLa_siArl8b_ratio <- subset(ratio_file, identifier == treatName)

# orga_intensityOnDetection_backsub
with(ratio_file, shapiro.test(intensity_ratio[identifier == controlName]))
with(ratio_file, shapiro.test(intensity_ratio[identifier == treatName]))

test <- wilcox.test(HeLa_scr_ratio$intensity_ratio, 
                    HeLa_siArl8b_ratio$intensity_ratio)
test$p.value
