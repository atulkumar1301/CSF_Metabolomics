#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)
library(scales)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
TABLE_Ab<-as.data.frame(matrix(ncol=10, nrow=376))
TABLE_Tau<-as.data.frame(matrix(ncol=10, nrow=376))
TABLE_Asyn<-as.data.frame(matrix(ncol=10, nrow=376))
TABLE_WMH<-as.data.frame(matrix(ncol=10, nrow=376))
names(TABLE_Ab)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC")
names(TABLE_Tau)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC")
names(TABLE_Asyn)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC")
names(TABLE_WMH)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")
d_1 <- rescale(df$samseg_wmhs_WMH_total_mm3, to = c (0, 1))
df <- add_column(df, DV = d_1, .after = 2)
df_2 <- df [,1:29]
j <- 1
for (i in colnames (df)) {
  if (i %in% colnames (df_2)) next
  N_P <- df[[i]]
  modeldata <- glm (N_P ~ 1, family=gaussian, data = df)
  model <- glm (N_P ~ Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + SAA_Status + DV + Age + Gender + Recruitment_Bias + mean_standardized_metabolomic_level, data = df, family=gaussian)
  lreg.or <-exp(cbind(OR = coef(model)))
  l0 <- deviance(modeldata);df0 <- df.residual(modeldata)
  l1 <- deviance(model);df1 <- df.residual(model)
  TABLE_Ab[j, 1] <- i
  TABLE_Ab[j,2] <- summary(model)$coefficients[2, "Estimate"]
  TABLE_Ab[j,3] <- lreg.or[2] # for odds ratio
  TABLE_Ab[j,4] <- summary(model)$coefficients[2, "Std. Error"]
  TABLE_Ab[j,5] <- summary(model)$coefficients[2, "Pr(>|t|)"]
  TABLE_Ab[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_Ab[j,7] <- confint(model) [2,1] 
  TABLE_Ab[j,8] <- confint(model) [2,2]
  TABLE_Ab[j,9] <- AIC (model)
  TABLE_Ab[j,10] <- BIC (model)
  
  TABLE_Tau[j, 1] <- i
  TABLE_Tau[j,2] <- summary(model)$coefficients[3, "Estimate"]
  TABLE_Tau[j,3] <- lreg.or[3] # for odds ratio
  TABLE_Tau[j,4] <- summary(model)$coefficients[3, "Std. Error"]
  TABLE_Tau[j,5] <- summary(model)$coefficients[3, "Pr(>|t|)"]
  TABLE_Tau[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_Tau[j,7] <- confint(model) [2,1] 
  TABLE_Tau[j,8] <- confint(model) [2,2]
  TABLE_Tau[j,9] <- AIC (model)
  TABLE_Tau[j,10] <- BIC (model)
  
  TABLE_Asyn[j, 1] <- i
  TABLE_Asyn[j,2] <- summary(model)$coefficients[4, "Estimate"]
  TABLE_Asyn[j,3] <- lreg.or[4] # for odds ratio
  TABLE_Asyn[j,4] <- summary(model)$coefficients[4, "Std. Error"]
  TABLE_Asyn[j,5] <- summary(model)$coefficients[4, "Pr(>|t|)"]
  TABLE_Asyn[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_Asyn[j,7] <- confint(model) [2,1] 
  TABLE_Asyn[j,8] <- confint(model) [2,2]
  TABLE_Asyn[j,9] <- AIC (model)
  TABLE_Asyn[j,10] <- BIC (model)
  
  TABLE_WMH[j, 1] <- i
  TABLE_WMH[j,2] <- summary(model)$coefficients[5, "Estimate"]
  TABLE_WMH[j,3] <- lreg.or[5] # for odds ratio
  TABLE_WMH[j,4] <- summary(model)$coefficients[5, "Std. Error"]
  TABLE_WMH[j,5] <- summary(model)$coefficients[5, "Pr(>|t|)"]
  TABLE_WMH[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_WMH[j,7] <- confint(model) [2,1] 
  TABLE_WMH[j,8] <- confint(model) [2,2]
  TABLE_WMH[j,9] <- AIC (model)
  TABLE_WMH[j,10] <- BIC (model)
  j <- j + 1
}
TABLE_Ab$P_Bonferroni <- p.adjust(TABLE_Ab$P, method = "bonferroni", n = length(TABLE_Ab$P))
TABLE_Ab$P_FDR <- p.adjust(TABLE_Ab$P, method = "fdr", n = length(TABLE_Ab$P))
write.table (TABLE_Ab, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/6_Result_Data_Analysis_Ab.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

TABLE_Tau$P_Bonferroni <- p.adjust(TABLE_Tau$P, method = "bonferroni", n = length(TABLE_Tau$P))
TABLE_Tau$P_FDR <- p.adjust(TABLE_Tau$P, method = "fdr", n = length(TABLE_Tau$P))
write.table (TABLE_Tau, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/6_Result_Data_Analysis_Taupet.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

TABLE_Asyn$P_Bonferroni <- p.adjust(TABLE_Asyn$P, method = "bonferroni", n = length(TABLE_Asyn$P))
TABLE_Asyn$P_FDR <- p.adjust(TABLE_Asyn$P, method = "fdr", n = length(TABLE_Asyn$P))
write.table (TABLE_Asyn, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/6_Result_Data_Analysis_Asyn.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

TABLE_WMH$P_Bonferroni <- p.adjust(TABLE_WMH$P, method = "bonferroni", n = length(TABLE_WMH$P))
TABLE_WMH$P_FDR <- p.adjust(TABLE_WMH$P, method = "fdr", n = length(TABLE_WMH$P))
write.table (TABLE_WMH, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/6_Result_Data_Analysis_WML.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
