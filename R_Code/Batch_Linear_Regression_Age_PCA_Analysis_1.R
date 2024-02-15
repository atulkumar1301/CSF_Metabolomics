#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)
library(scales)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=10, nrow=369))
names(TABLE)<-c("Biomarker", "Effect", "OR","SE", "P", "R2", "L95", "U95", "AIC", "BIC")
df <- fread (file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/11_Data_Analysis_1_PCA_X.txt")
df_2 <- df [,1:29]
j <- 1
for (i in colnames (df)) {
  if (i %in% colnames (df_2)) next
  modeldata <- glm (Age ~ 1, family=gaussian, data = df)
  N_P <- df[[i]]
  model <- glm (Age ~ N_P, data = df, family=gaussian)
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
  j <- j + 1
}
write.table (TABLE, (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/12_Result_Data_Analysis_Age_Raw.txt")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
