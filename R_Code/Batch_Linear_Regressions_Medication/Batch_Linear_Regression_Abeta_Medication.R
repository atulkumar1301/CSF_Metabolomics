
#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)
library(scales)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=11, nrow=376))
names(TABLE)<-c("Protein", "Effect", "OR","SE", "P", "R2", "t Value", "L95", "U95", "AIC", "BIC")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")
df_2 <- df [,1:29]
j <- 1
for (i in colnames (df)) {
  if (i %in% colnames (df_2)) next
  N_P <- df[[i]]
  modeldata <- glm (N_P ~ 1, family=gaussian, data = df)
  model <- glm (N_P ~ Abnormal_CSF_Ab42_Ab40_Ratio + Age + Gender + Recruitment_Bias + mean_standardized_metabolomic_level + AD_PD_Medication, data = df, family=gaussian)
  lreg.or <-exp(cbind(OR = coef(model)))
  TABLE[j, 1] <- i
  TABLE[j,2] <- summary(model)$coefficients[2, "Estimate"]
  TABLE[j,3] <- lreg.or[2]
  TABLE[j,4] <- summary(model)$coefficients[2, "Std. Error"]
  TABLE[j,5] <- summary(model)$coefficients[2, "Pr(>|t|)"]
  l0 <- deviance(modeldata);df0 <- df.residual(modeldata)
  l1 <- deviance(model);df1 <- df.residual(model)
  TABLE[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  obs = model$y
  pred = model$fitted.values
  TABLE[j,7] <- summary(model)$coefficients[2, "t value"]
  TABLE[j,8] <- confint(model) [2,1] 
  TABLE[j,9] <- confint(model) [2,2]
  TABLE[j,10] <- AIC (model)
  TABLE[j,11] <- BIC (model)
  j <- j + 1
}
TABLE$P_Bonferroni <- p.adjust(TABLE$P, method = "bonferroni", n = length(TABLE$P))
TABLE$P_FDR <- p.adjust(TABLE$P, method = "fdr", n = length(TABLE$P))
write.table (TABLE, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/A_Beta/Medication/6_Result_Data_Analysis_A_beta_Medication.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
