#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)
library(scales)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=11, nrow=71))
names(TABLE)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Full_Data_Zscore_Average.txt")
j <- 1
d_1 <- rescale(df$samseg_wmhs_WMH_total_mm3, to = c (0, 1))
#d_2 <- rescale(df$icv_mm3, to = c (0, 1))
df <- add_column(df, DV = d_1, .after = 2)
#df <- add_column(df, Cov = d_2, .after = 2)
df_2 <- df [,1:29]
for (i in colnames (df)) {
  if (i %in% colnames (df_2)) next
  N_P <- df[[i]]
  modeldata <- glm (N_P ~ 1, family=gaussian, data = df)
  model <- glm (N_P ~ DV + Age + Gender + Recruitment_Bias + icv_mm3 + mean_standardized_metabolomic_level, data = df, family=gaussian)
  lreg.or <-exp(cbind(OR = coef(model)))
  TABLE[j, 1] <- i
  TABLE[j,2] <- summary(model)$coefficients[2, "Estimate"]
  TABLE[j,3] <- lreg.or[2] # for odds ratio
  TABLE[j,4] <- summary(model)$coefficients[2, "Std. Error"]
  TABLE[j,5] <- summary(model)$coefficients[2, "Pr(>|t|)"]
  l0 <- deviance(modeldata);df0 <- df.residual(modeldata)
  l1 <- deviance(model);df1 <- df.residual(model)
  TABLE[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE[j,7] <- confint(model) [2,1] 
  TABLE[j,8] <- confint(model) [2,2]
  TABLE[j,9] <- AIC (model)
  TABLE[j,10] <- BIC (model)
  TABLE[j,11] <- summary(model)$coefficients[2, "t value"]
  j <- j + 1
}
TABLE$P_Bonferroni <- p.adjust(TABLE$P, method = "bonferroni", n = length(TABLE$P))
TABLE$P_FDR <- p.adjust(TABLE$P, method = "fdr", n = length(TABLE$P))
write.table (TABLE, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/WML/6_Result_Data_Analysis_WML.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
