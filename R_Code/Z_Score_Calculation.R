#! /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/bin/Rscript
library (data.table)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
df <- fread (file = paste0("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Data_Sub_Pathway/", args [1]))
df1 <- df [,-1]
df2 <- as.data.frame (sapply(df1, function(df1) (df1-mean(df1))/sd(df1)))
df2$Average <- rowMeans(df2)
write.table (df2, file = paste0("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Zscore_Data/", args[1]), sep = "\t", row.names = FALSE, quote=FALSE)
