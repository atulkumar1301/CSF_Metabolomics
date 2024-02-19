library(ggVennDiagram)
library(data.table)
library(ggplot2)
library(tidyverse)
### Abeta Medication, Abeta
df_abeta <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/A_Beta/8_Full_Result_Data_Analysis_A_beta.txt")
df_abeta_Sig_Pos <- df_abeta[df_abeta$P_FDR <= 0.05 & df_abeta$Effect > 0]
df_abeta_Sig_Neg <- df_abeta[df_abeta$P_FDR <= 0.05 & df_abeta$Effect < 0]
abeta_Sig_Pos <- df_abeta_Sig_Pos$Protein
abeta_Sig_Neg <- df_abeta_Sig_Neg$Protein

df_abeta_Medication <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/A_Beta/Medication/8_Full_Result_Data_Analysis_A_beta_Medication.txt")
df_abeta_Medication_Sig_Pos <- df_abeta_Medication[df_abeta_Medication$P_FDR <= 0.05 & df_abeta_Medication$Effect > 0]
df_abeta_Medication_Sig_Neg <- df_abeta_Medication[df_abeta_Medication$P_FDR <= 0.05 & df_abeta_Medication$Effect < 0]
abeta_Medication_Sig_Pos <- df_abeta_Medication_Sig_Pos$Protein
abeta_Medication_Sig_Neg <- df_abeta_Medication_Sig_Neg$Protein

## Positive
venn_Pos <- list (abeta_Sig_Pos, abeta_Medication_Sig_Pos)
p <- ggVennDiagram (venn_Pos, label_alpha = 0,category.names = c("Abeta", "Medication"))
p <- p + ggplot2::scale_fill_gradient(low="#D55E00",high = "#009E73") +
  ggtitle (expression ("Increased level of Metabolites in A"*beta*" Status vs A"*beta*" Status adjusted for Medication"))
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
venn_Neg <- list (abeta_Sig_Neg, abeta_Medication_Sig_Neg)
p1 <- ggVennDiagram (venn_Neg, label_alpha = 0,category.names = c("Abeta", "Medication"))
p1 <- p1 + ggplot2::scale_fill_gradient(low="#E69F00",high = "#56B4E9") +
  ggtitle (expression ("Decreased level of Metabolites in A"*beta*" Status vs A"*beta*" Status adjusted for Medication"))
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
