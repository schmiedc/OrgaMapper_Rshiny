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