library(ggVennDiagram)
library(data.table)
library(ggplot2)
library(tidyverse)
### Tau Medication, Tau
df_Tau <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/TauPET/8_Full_Result_Data_Analysis_Taupet.txt")
df_Tau_Sig_Pos <- df_Tau[df_Tau$P_FDR <= 0.05 & df_Tau$Effect > 0]
df_Tau_Sig_Neg <- df_Tau[df_Tau$P_FDR <= 0.05 & df_Tau$Effect < 0]
Tau_Sig_Pos <- df_Tau_Sig_Pos$Biomarker
Tau_Sig_Neg <- df_Tau_Sig_Neg$Biomarker

df_Tau_Medication <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/TauPET//Medication/8_Full_Result_Data_Analysis_Taupet_Medication.txt")
df_Tau_Medication_Sig_Pos <- df_Tau_Medication[df_Tau_Medication$P_FDR <= 0.05 & df_Tau_Medication$Effect > 0]
df_Tau_Medication_Sig_Neg <- df_Tau_Medication[df_Tau_Medication$P_FDR <= 0.05 & df_Tau_Medication$Effect < 0]
Tau_Medication_Sig_Pos <- df_Tau_Medication_Sig_Pos$Biomarker
Tau_Medication_Sig_Neg <- df_Tau_Medication_Sig_Neg$Biomarker

df_Tau_Ab40 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/TauPET/Ab40/7_Full_Result_Data_Analysis_Taupet_Ab40.txt")
df_Tau_Ab40_Sig_Pos <- df_Tau_Ab40[df_Tau_Ab40$P_FDR <= 0.05 & df_Tau_Ab40$Effect > 0]
df_Tau_Ab40_Sig_Neg <- df_Tau_Ab40[df_Tau_Ab40$P_FDR <= 0.05 & df_Tau_Ab40$Effect < 0]
Tau_Ab40_Sig_Pos <- df_Tau_Ab40_Sig_Pos$Biomarker
Tau_Ab40_Sig_Neg <- df_Tau_Ab40_Sig_Neg$Biomarker

## Positive
venn_Pos <- list (Tau_Sig_Pos, Tau_Medication_Sig_Pos, Tau_Ab40_Sig_Pos)
p <- ggVennDiagram (venn_Pos, label_alpha = 0,category.names = c("TauPET", "Medication", "Ab40"))
p <- p + ggplot2::scale_fill_gradient(low="#D55E00",high = "#009E73") +
  ggtitle (expression ("Increased level of Metabolites in TauPET, TauPET adjusted for Medication and TauPET adjusted for A"*beta*"40"))
p <- p + theme(
  plot.title = element_text(family = "serif", size=18, face = "bold"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  legend.title = element_text(family = "serif", size=16),
  legend.text = element_text(family = "serif", size=16),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank())
p

## Negative
venn_Neg <- list (Tau_Sig_Neg, Tau_Medication_Sig_Neg, Tau_Ab40_Sig_Neg)
p1 <- ggVennDiagram (venn_Neg, label_alpha = 0,category.names = c("TauPET", "Medication", "Ab40"))
p1 <- p1 + ggplot2::scale_fill_gradient(low="#E69F00",high = "#56B4E9") +
  ggtitle (expression ("Decreased level of Metabolites in TauPET, TauPET adjusted for Medication and TauPET adjusted for A"*beta*"40"))
p1 <- p1 + theme(
  plot.title = element_text(family = "serif", size=18, face = "bold"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  legend.title = element_text(family = "serif", size=16),
  legend.text = element_text(family = "serif", size=16),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank())
p1
