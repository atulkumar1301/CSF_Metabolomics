#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(scales)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=12, nrow=376))
names(TABLE)<-c("Model", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC", "t Value")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt")

j = 1
modeldata_Ab <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ 1, family=binomial (link = 'logit'), data = df)
model_Ab <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ mean_standardized_metabolomic_level + Age + Gender + Recruitment_Bias, data = df, family=binomial (link = 'logit'))
TABLE[j, 1] <- "Abeta"
TABLE[j, 2] <- summary(model_Ab)$coefficients[2, "Estimate"]
#TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Ab)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Ab)$coefficients[2, "Pr(>|z|)"]
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
modeldata_Ab_Ad <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ 1, family=binomial (link = 'logit'), data = df)
model_Ab_Ad <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ mean_standardized_metabolomic_level + Age + Gender + Recruitment_Bias, data = df, family=binomial (link = 'logit'))
TABLE[j, 1] <- "Abeta Adjusted"
TABLE[j, 2] <- summary(model_Ab_Ad)$coefficients[2, "Estimate"]
#TABLE[j, 3] <- lreg.or[2]
TABLE[j, 4] <- summary(model_Ab_Ad)$coefficients[2, "Std. Error"]
TABLE[j, 5] <- summary(model_Ab_Ad)$coefficients[2, "Pr(>|z|)"]
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




TABLE$P_Bonferroni <- p.adjust(TABLE$P, method = "bonferroni", n = length(TABLE$P))
TABLE$P_FDR <- p.adjust(TABLE$P, method = "fdr", n = length(TABLE$P))
write.table (TABLE, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/A_Beta/6_Result_Data_Analysis_A_beta.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
