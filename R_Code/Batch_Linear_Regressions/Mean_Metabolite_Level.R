#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(scales)
library(tibble)
library(ROCR)
library(pROC)
library(RNOmni)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=11, nrow=8))
names(TABLE)<-c("Model", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")

j = 1
modeldata_Ab <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Ab <- glm (mean_standardized_metabolomic_level ~ Abnormal_CSF_Ab42_Ab40_Ratio + Age + Gender + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Ab)))
TABLE[j, 1] <- "Abeta"
TABLE[j, 2] <- summary(model_Ab)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Ab)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Ab)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Ab);df0 <- df.residual(modeldata_Ab)
l1 <- deviance(model_Ab);df1 <- df.residual(model_Ab)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Ab$y
pred = model_Ab$fitted.values
TABLE[j, 7] <- confint(model_Ab) [2,1] 
TABLE[j, 8] <- confint(model_Ab) [2,2]
TABLE[j, 9] <- AIC (model_Ab)
TABLE[j, 10] <- BIC (model_Ab)
TABLE[j, 11] <- summary(model_Ab)$coefficients[2, "t value"]

j = 2
modeldata_Ab_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Ab_Ad <- glm (mean_standardized_metabolomic_level ~ Abnormal_CSF_Ab42_Ab40_Ratio + Age + Gender + Recruitment_Bias + SAA_Status + samseg_wmhs_WMH_total_mm3 + tnic_cho_com_I_IV, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Ab_Ad)))
TABLE[j, 1] <- "Abeta Adjusted"
TABLE[j, 2] <- summary(model_Ab_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Ab_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Ab_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Ab_Ad);df0 <- df.residual(modeldata_Ab_Ad)
l1 <- deviance(model_Ab_Ad);df1 <- df.residual(model_Ab_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Ab_Ad$y
pred = model_Ab_Ad$fitted.values
TABLE[j, 7] <- confint(model_Ab_Ad) [2,1] 
TABLE[j, 8] <- confint(model_Ab_Ad) [2,2]
TABLE[j, 9] <- AIC (model_Ab_Ad)
TABLE[j, 10] <- BIC (model_Ab_Ad)
TABLE[j, 11] <- summary(model_Ab_Ad)$coefficients[2, "t value"]

j = 3
modeldata_Tau <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Tau <- glm (mean_standardized_metabolomic_level ~ tnic_cho_com_I_IV + Age + Gender + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Tau)))
TABLE[j, 1] <- "Tau PET"
TABLE[j, 2] <- summary(model_Tau)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Tau)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Tau)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Tau);df0 <- df.residual(modeldata_Tau)
l1 <- deviance(model_Tau);df1 <- df.residual(model_Tau)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Tau$y
pred = model_Tau$fitted.values
TABLE[j, 7] <- confint(model_Tau) [2,1] 
TABLE[j, 8] <- confint(model_Tau) [2,2]
TABLE[j, 9] <- AIC (model_Tau)
TABLE[j, 10] <- BIC (model_Tau)
TABLE[j, 11] <- summary(model_Tau)$coefficients[2, "t value"]

j = 4
modeldata_Tau_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Tau_Ad <- glm (mean_standardized_metabolomic_level ~ tnic_cho_com_I_IV + Age + Gender + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + SAA_Status + samseg_wmhs_WMH_total_mm3, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Tau_Ad)))
TABLE[j, 1] <- "Tau PET Adjusted"
TABLE[j, 2] <- summary(model_Tau_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Tau_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Tau_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Tau_Ad);df0 <- df.residual(modeldata_Tau_Ad)
l1 <- deviance(model_Tau_Ad);df1 <- df.residual(model_Tau_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Tau_Ad$y
pred = model_Tau_Ad$fitted.values
TABLE[j, 7] <- confint(model_Tau_Ad) [2,1] 
TABLE[j, 8] <- confint(model_Tau_Ad) [2,2]
TABLE[j, 9] <- AIC (model_Tau_Ad)
TABLE[j, 10] <- BIC (model_Tau_Ad)
TABLE[j, 11] <- summary(model_Tau_Ad)$coefficients[2, "t value"]

j = 5
modeldata_SAA <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_SAA <- glm (mean_standardized_metabolomic_level ~ SAA_Status + Age + Gender + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_SAA)))
TABLE[j, 1] <- "SAA"
TABLE[j, 2] <- summary(model_SAA)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_SAA)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_SAA)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_SAA);df0 <- df.residual(modeldata_SAA)
l1 <- deviance(model_SAA);df1 <- df.residual(model_SAA)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_SAA$y
pred = model_SAA$fitted.values
TABLE[j, 7] <- confint(model_SAA) [2,1] 
TABLE[j, 8] <- confint(model_SAA) [2,2]
TABLE[j, 9] <- AIC (model_SAA)
TABLE[j, 10] <- BIC (model_SAA)
TABLE[j, 11] <- summary(model_SAA)$coefficients[2, "t value"]

j = 6
modeldata_SAA_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_SAA_Ad <- glm (mean_standardized_metabolomic_level ~ SAA_Status + Age + Gender + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_SAA_Ad)))
TABLE[j, 1] <- "SAA Adjusted"
TABLE[j, 2] <- summary(model_SAA_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_SAA_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_SAA_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_SAA_Ad);df0 <- df.residual(modeldata_SAA_Ad)
l1 <- deviance(model_SAA_Ad);df1 <- df.residual(model_SAA_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_SAA_Ad$y
pred = model_SAA_Ad$fitted.values
TABLE[j, 7] <- confint(model_SAA_Ad) [2,1] 
TABLE[j, 8] <- confint(model_SAA_Ad) [2,2]
TABLE[j, 9] <- AIC (model_SAA_Ad)
TABLE[j, 10] <- BIC (model_SAA_Ad)
TABLE[j, 11] <- summary(model_SAA_Ad)$coefficients[2, "t value"]

j = 7
d_1 <- rescale(df$samseg_wmhs_WMH_total_mm3, to = c (0, 1))
modeldata_WML <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_WML <- glm (mean_standardized_metabolomic_level ~ d_1 + Age + Gender + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_WML)))
TABLE[j, 1] <- "WML"
TABLE[j, 2] <- summary(model_WML)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_WML)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_WML)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_WML);df0 <- df.residual(modeldata_WML)
l1 <- deviance(model_WML);df1 <- df.residual(model_WML)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_WML$y
pred = model_WML$fitted.values
TABLE[j, 7] <- confint(model_WML) [2,1] 
TABLE[j, 8] <- confint(model_WML) [2,2]
TABLE[j, 9] <- AIC (model_WML)
TABLE[j, 10] <- BIC (model_WML)
TABLE[j, 11] <- summary(model_WML)$coefficients[2, "t value"]

j = 8
modeldata_WML_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_WML_Ad <- glm (mean_standardized_metabolomic_level ~ d_1 + Age + Gender + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + SAA_Status, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_WML_Ad)))
TABLE[j, 1] <- "WML Adjusted"
TABLE[j, 2] <- summary(model_WML_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_WML_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_WML_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_WML_Ad);df0 <- df.residual(modeldata_WML_Ad)
l1 <- deviance(model_WML_Ad);df1 <- df.residual(model_WML_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_WML_Ad$y
pred = model_WML_Ad$fitted.values
TABLE[j, 7] <- confint(model_WML_Ad) [2,1] 
TABLE[j, 8] <- confint(model_WML_Ad) [2,2]
TABLE[j, 9] <- AIC (model_WML_Ad)
TABLE[j, 10] <- BIC (model_WML_Ad)
TABLE[j, 11] <- summary(model_WML_Ad)$coefficients[2, "t value"]

write.table (TABLE, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Mean_Metabolite_Level_Pathology/Result_Mean_Metabolite_Level_Pathology.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
