library(data.table)
library(ggpubr)
library(ggplot2)
library(ggcorrplot)
library(corrr)
library(dplyr)
##########
df_Metabolites <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")
Mat_Meta <- df_Metabolites [,-1:-34]
Mat_Meta_1 <- as.matrix(Mat_Meta)
Meta_Mat <- correlate (Mat_Meta_1)
Meta_Mat <- Meta_Mat [,-1]
Meta_Mat_p <- ggcorrplot(Meta_Mat, outline.col = "white", lab = FALSE, hc.order = FALSE) + ggtitle ("Correlation Matrix Metabolites")
Meta_Mat_p <- Meta_Mat_p + scale_x_continuous(labels = c ())
Meta_Mat_p  <- Meta_Mat_p + scale_y_discrete(labels = c ())
Meta_Mat_p <- Meta_Mat_p + theme(
    plot.title = element_text(family = "serif", size=18, face = "bold"),
    #axis.title.x = element_text(family = "serif", size=16),
    #axis.title.y = element_text(family = "serif", size=16),
    axis.text.x = element_text(family = "serif", size=12),
    axis.text.y = element_text(family = "serif", size=12),
    legend.title = element_text(family = "serif", size=12),
    legend.text = element_text(family = "serif", size=12),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())
