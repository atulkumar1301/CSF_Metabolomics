#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)

library(scales)
TABLE_Age <- as.data.frame(matrix(ncol=11, nrow=71))
names(TABLE_Age) <- c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")

TABLE_Gender <-as.data.frame(matrix(ncol=11, nrow=71))
names(TABLE_Gender) <- c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")

TABLE_CU_CI <- as.data.frame(matrix(ncol=11, nrow=71))
names(TABLE_CU_CI) <- c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")

df <- fread (file = "~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Full_Data_Zscore_Average.txt")
df_2 <- df [,1:28]
j <- 1
for (i in colnames (df)) {
  if (i %in% colnames (df_2)) next
  N_P <- df[[i]]
  modeldata <- glm (N_P ~ 1, family=gaussian, data = df)
  model <- glm (N_P ~ Age + Gender + Recruitment_Bias + mean_standardized_metabolomic_level, data = df, family=gaussian)
  lreg.or <-exp(cbind(OR = coef(model)))
  l0 <- deviance(modeldata);df0 <- df.residual(modeldata)
  l1 <- deviance(model);df1 <- df.residual(model)
  
  TABLE_Age[j, 1] <- i
  TABLE_Age[j,2] <- summary(model)$coefficients[2, "Estimate"]
  TABLE_Age[j,3] <- lreg.or[2] # for odds ratio
  TABLE_Age[j,4] <- summary(model)$coefficients[2, "Std. Error"]
  TABLE_Age[j,5] <- summary(model)$coefficients[2, "Pr(>|t|)"]
  TABLE_Age[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_Age[j,7] <- confint(model) [2,1] 
  TABLE_Age[j,8] <- confint(model) [2,2]
  TABLE_Age[j,9] <- AIC (model)
  TABLE_Age[j,10] <- BIC (model)
  TABLE_Age[j,11] <- summary(model)$coefficients[2, "t value"]
  
  TABLE_Gender[j, 1] <- i
  TABLE_Gender[j,2] <- summary(model)$coefficients[3, "Estimate"]
  TABLE_Gender[j,3] <- lreg.or[2] # for odds ratio
  TABLE_Gender[j,4] <- summary(model)$coefficients[3, "Std. Error"]
  TABLE_Gender[j,5] <- summary(model)$coefficients[3, "Pr(>|t|)"]
  TABLE_Gender[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_Gender[j,7] <- confint(model) [3,1] 
  TABLE_Gender[j,8] <- confint(model) [3,2]
  TABLE_Gender[j,9] <- AIC (model)
  TABLE_Gender[j,10] <- BIC (model)
  TABLE_Gender[j,11] <- summary(model)$coefficients[3, "t value"]
  
  TABLE_CU_CI[j, 1] <- i
  TABLE_CU_CI[j,2] <- summary(model)$coefficients[4, "Estimate"]
  TABLE_CU_CI[j,3] <- lreg.or[2] # for odds ratio
  TABLE_CU_CI[j,4] <- summary(model)$coefficients[4, "Std. Error"]
  TABLE_CU_CI[j,5] <- summary(model)$coefficients[4, "Pr(>|t|)"]
  TABLE_CU_CI[j,6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
  TABLE_CU_CI[j,7] <- confint(model) [4,1] 
  TABLE_CU_CI[j,8] <- confint(model) [4,2]
  TABLE_CU_CI[j,9] <- AIC (model)
  TABLE_CU_CI[j,10] <- BIC (model)
  TABLE_CU_CI[j,11] <- summary(model)$coefficients[4, "t value"]
  
  j <- j + 1
}
TABLE_Age$P_Bonferroni <- p.adjust(TABLE_Age$P, method = "bonferroni", n = length(TABLE_Age$P))
TABLE_Age$P_FDR <- p.adjust(TABLE_Age$P, method = "fdr", n = length(TABLE_Age$P))
write.table (TABLE_Age, (file = paste0 ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/All_Combined/6_Result_Data_Analysis_Age.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

TABLE_Gender$P_Bonferroni <- p.adjust(TABLE_Gender$P, method = "bonferroni", n = length(TABLE_Gender$P))
TABLE_Gender$P_FDR <- p.adjust(TABLE_Gender$P, method = "fdr", n = length(TABLE_Gender$P))
write.table (TABLE_Gender, (file = paste0 ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/6_Result_Data_Analysis_Gender.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

TABLE_CU_CI$P_Bonferroni <- p.adjust(TABLE_CU_CI$P, method = "bonferroni", n = length(TABLE_CU_CI$P))
TABLE_CU_CI$P_FDR <- p.adjust(TABLE_CU_CI$P, method = "fdr", n = length(TABLE_CU_CI$P))
write.table (TABLE_CU_CI, (file = paste0 ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/All_Combined/6_Result_Data_Analysis_CU_CI.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
