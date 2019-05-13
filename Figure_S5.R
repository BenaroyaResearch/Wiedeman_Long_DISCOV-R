

library(tidyverse)
library(readxl)
library(hexbin)

specificity_phenos = read_excel("allSubjectsSummary_Compiled_10K_Samples1_EJ_Flow_Comparison_181016.xlsx")

# Jittered dot plot
ggplot(specificity_phenos, aes(x = `Flow specificities`, y = `CyTOF phenotypes`)) +
  geom_point(position=position_jitter(h=0.2, w=0.2), alpha = 0.7, size = 2) +
  stat_smooth(method = "lm", color = "black", se = F) + ylim(0,3.5)

# 2D histogram (hex)
ggplot(specificity_phenos, aes(x = `Flow specificities`, y = `CyTOF phenotypes`)) +
  stat_binhex(binwidth = c(1,1)) + 
  stat_smooth(method = "lm", color = "deeppink4", se = F, fullrange = T, size = 1.25) + 
  ylim(-0.5,3.75) + xlim(-0.5,3.75)
  
# Contingency table
table(specificity_phenos$`CyTOF phenotypes`, specificity_phenos$`Flow specificities`)
