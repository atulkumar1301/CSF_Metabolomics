#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(scales)
library(tibble)
library(ROCR)
library(pROC)
library(RNOmni)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=11, nrow=14))
names(TABLE)<-c("Model", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")
df <- fread (file = "~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")

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

j = 9
modeldata_Age <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Age <- glm (mean_standardized_metabolomic_level ~ Age + Gender + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Age)))
TABLE[j, 1] <- "Age"
TABLE[j, 2] <- summary(model_Age)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Age)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Age)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Age);df0 <- df.residual(modeldata_Age)
l1 <- deviance(model_Age);df1 <- df.residual(model_Age)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Age$y
pred = model_Age$fitted.values
TABLE[j, 7] <- confint(model_Age) [2,1] 
TABLE[j, 8] <- confint(model_Age) [2,2]
TABLE[j, 9] <- AIC (model_Age)
TABLE[j, 10] <- BIC (model_Age)
TABLE[j, 11] <- summary(model_Age)$coefficients[2, "t value"]

j = 10
modeldata_Age_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Age_Ad <- glm (mean_standardized_metabolomic_level ~ Age + Gender + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + SAA_Status, samseg_wmhs_WMH_total_mm3, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Age_Ad)))
TABLE[j, 1] <- "Age Adjusted"
TABLE[j, 2] <- summary(model_Age_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Age_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Age_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Age_Ad);df0 <- df.residual(modeldata_Age_Ad)
l1 <- deviance(model_Age_Ad);df1 <- df.residual(model_Age_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Age_Ad$y
pred = model_Age_Ad$fitted.values
TABLE[j, 7] <- confint(model_Age_Ad) [2,1] 
TABLE[j, 8] <- confint(model_Age_Ad) [2,2]
TABLE[j, 9] <- AIC (model_Age_Ad)
TABLE[j, 10] <- BIC (model_Age_Ad)
TABLE[j, 11] <- summary(model_Age_Ad)$coefficients[2, "t value"]

j = 11
modeldata_Gender <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Gender <- glm (mean_standardized_metabolomic_level ~ Gender + Age + Recruitment_Bias, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Gender)))
TABLE[j, 1] <- "Gender"
TABLE[j, 2] <- summary(model_Gender)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Gender)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Gender)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Gender);df0 <- df.residual(modeldata_Gender)
l1 <- deviance(model_Gender);df1 <- df.residual(model_Gender)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Gender$y
pred = model_Gender$fitted.values
TABLE[j, 7] <- confint(model_Gender) [2,1] 
TABLE[j, 8] <- confint(model_Gender) [2,2]
TABLE[j, 9] <- AIC (model_Gender)
TABLE[j, 10] <- BIC (model_Gender)
TABLE[j, 11] <- summary(model_Gender)$coefficients[2, "t value"]

j = 12
modeldata_Gender_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_Gender_Ad <- glm (mean_standardized_metabolomic_level ~ Gender + Age + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + SAA_Status, samseg_wmhs_WMH_total_mm3, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_Gender_Ad)))
TABLE[j, 1] <- "Gender Adjusted"
TABLE[j, 2] <- summary(model_Gender_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Gender_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Gender_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_Gender_Ad);df0 <- df.residual(modeldata_Gender_Ad)
l1 <- deviance(model_Gender_Ad);df1 <- df.residual(model_Gender_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_Gender_Ad$y
pred = model_Gender_Ad$fitted.values
TABLE[j, 7] <- confint(model_Gender_Ad) [2,1] 
TABLE[j, 8] <- confint(model_Gender_Ad) [2,2]
TABLE[j, 9] <- AIC (model_Gender_Ad)
TABLE[j, 10] <- BIC (model_Gender_Ad)
TABLE[j, 11] <- summary(model_Gender_Ad)$coefficients[2, "t value"]

j = 13
modeldata_CU_CI <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_CU_CI <- glm (mean_standardized_metabolomic_level ~ Recruitment_Bias + Gender + Age, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_CU_CI)))
TABLE[j, 1] <- "CU_CI"
TABLE[j, 2] <- summary(model_CU_CI)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_CU_CI)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_CU_CI)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_CU_CI);df0 <- df.residual(modeldata_CU_CI)
l1 <- deviance(model_CU_CI);df1 <- df.residual(model_CU_CI)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_CU_CI$y
pred = model_CU_CI$fitted.values
TABLE[j, 7] <- confint(model_CU_CI) [2,1] 
TABLE[j, 8] <- confint(model_CU_CI) [2,2]
TABLE[j, 9] <- AIC (model_CU_CI)
TABLE[j, 10] <- BIC (model_CU_CI)
TABLE[j, 11] <- summary(model_CU_CI)$coefficients[2, "t value"]

j = 14
modeldata_CU_CI_Ad <- glm (mean_standardized_metabolomic_level ~ 1, family=gaussian, data = df)
model_CU_CI_Ad <- glm (mean_standardized_metabolomic_level ~ Recruitment_Bias + Gender + Age + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + SAA_Status, samseg_wmhs_WMH_total_mm3, data = df, family=gaussian)
lreg.or <-exp(cbind(OR = coef(model_CU_CI_Ad)))
TABLE[j, 1] <- "CU_CI Adjusted"
TABLE[j, 2] <- summary(model_CU_CI_Ad)$coefficients[2, "Estimate"]
TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_CU_CI_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_CU_CI_Ad)$coefficients[2, "Pr(>|t|)"]
l0 <- deviance(modeldata_CU_CI_Ad);df0 <- df.residual(modeldata_CU_CI_Ad)
l1 <- deviance(model_CU_CI_Ad);df1 <- df.residual(model_CU_CI_Ad)
TABLE[j, 6] <- (1 - exp((l1 - l0)/nrow(df)))/(1 - exp( - l0/nrow(df)))        #Nagelkerke
obs = model_CU_CI_Ad$y
pred = model_CU_CI_Ad$fitted.values
TABLE[j, 7] <- confint(model_CU_CI_Ad) [2,1] 
TABLE[j, 8] <- confint(model_CU_CI_Ad) [2,2]
TABLE[j, 9] <- AIC (model_CU_CI_Ad)
TABLE[j, 10] <- BIC (model_CU_CI_Ad)
TABLE[j, 11] <- summary(model_CU_CI_Ad)$coefficients[2, "t value"]

write.table (TABLE, (file = paste0 ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Mean_Metabolite_Level_Pathology/Result_Mean_Metabolite_Level_Pathology.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
