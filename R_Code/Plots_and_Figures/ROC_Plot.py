#Plotting ROC Curve
library(data.table)
library(ROCR)
library(pROC)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/PD_CRH/Data_SAA-Ctrl_vs_SAA+LBD.txt")
modeldata <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ 1, family=binomial (link = 'logit'), data = df)
model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `carboxyethyl-GABA` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")




#Box Plot
library(data.table)
par(mfrow=c(2,2)) # for multiple plots
df <- fread ("Full_Data_Metabolites.txt")
par(cex.axis = 1.5)
av <- aov(`carboxyethyl-GABA` ~ Abnormal_CSF_Ab42_Ab40_Ratio, data = df)
summary (av)
boxplot(`carboxyethyl-GABA` ~ Abnormal_CSF_Ab42_Ab40_Ratio, data = df, lwd = 2, ylab = 'Carboxyethyl-GABA', notch = TRUE, cex.lab = 1.5, main = "A) Carboxyethyl-GABA", cex.main = 1.5)
stripchart(`carboxyethyl-GABA` ~ Abnormal_CSF_Ab42_Ab40_Ratio, vertical = TRUE, data = df, 
           method = "jitter", add = TRUE, pch = 20, col = c('#E69F00', '#CC79A7'), cex = 2, jitter = 0.3)
mtext("Top six metabolites associted with Ab", side = 3, line = -1.75, outer = TRUE, cex = 1.75, col = '#D55E00')

df_1 <- fread ("/Volumes/ATUL_6TB/Work/Projects/Nicholas_Ashton/Data_SAA-Ctrl_vs_SAA+LBD_Denovo.txt")
par(cex.axis=1.5)
boxplot(CRH__Oncology_o2csf__NPX ~ Group, data = df_1, lwd = 2, ylab = 'CSF CRH', notch = TRUE, cex.lab = 1.5, main = "B) SAA- CTR vs De Novo SAA+ LBD", cex.main = 1.5)
stripchart(CRH__Oncology_o2csf__NPX ~ Group, vertical = TRUE, data = df_1, 
           method = "jitter", add = TRUE, pch = 20, col = c('#E69F00', '#CC79A7'), cex = 2, jitter = 0.3)

df_2 <- fread ("/Volumes/ATUL_6TB/Work/Projects/Nicholas_Ashton/Data_SAA-Ctrl_vs_SAA+Ctrl.txt")
par(cex.axis=1.5)
boxplot(CRH__Oncology_o2csf__NPX ~ Group, data = df_2, lwd = 2, ylab = 'CSF CRH', notch = TRUE, cex.lab = 1.5, main = "C) SAA- CTR vs SAA+ CTR", cex.main = 1.5)
stripchart(CRH__Oncology_o2csf__NPX ~ Group, vertical = TRUE, data = df_2, 
           method = "jitter", add = TRUE, pch = 20, col = c('#E69F00', '#CC79A7'), cex = 2, jitter = 0.3)


df_3 <- fread ("/Volumes/ATUL_6TB/Work/Projects/Nicholas_Ashton/Data_SAA-AD_FTD_VaD_vs_SAA+_LBD.txt")
par(cex.axis=1.5)
boxplot(CRH__Oncology_o2csf__NPX ~ Group, data = df_3, lwd = 2, ylab = 'CSF CRH', notch = TRUE, cex.lab = 1.5, main = "D) SAA- AD/FTD/VaD vs SAA+ LBD", cex.main = 1.5, col.main = "#0072B2")
stripchart(CRH__Oncology_o2csf__NPX ~ Group, vertical = TRUE, data = df_3, 
           method = "jitter", add = TRUE, pch = 20, col = c('#E69F00', '#0072B2'), cex = 2, jitter = 0.3)



#Plotting ROC Curve
library(data.table)
library(ROCR)
library(pROC)
library(RNOmni)
library(caret)
library(ggplot2)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/PD_CRH/Atypical_Analysis/Data_SAA-AD_FTD_VaD_vs_Atypical.txt")
modeldata <- glm (Diagnosis ~ 1, family=binomial (link = 'logit'), data = df)
N_P <- df$CRH__Oncology_o2csf__NPX
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_1 <- plot.roc(obs, pred,  col="#F0E442", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = FALSE, print.auc.x = 0.2, print.auc.y = 0.8, print.auc.pattern = "AUC = %.2f", print.auc.cex=1.5, cex.lab = 1.5, main = "ROC Curve: SAA- AD/FTD/VaD vs Atypical PS", cex.main = 1.5)
legend("bottomright", legend=c("CRH", "DDC", "NFL", "CRH/DDC Ratio", "CRH/NfL Ratio", "(CRH + NfL)/DDC Ratio", "CRH/(DDC + NfL) Ratio"), col=c("#F0E442", "#E69F00", "#009E73", "#CC79A7", "#0072B2", "#56B4E9", "#D55E00"), lwd=2, cex = 1.5)

N_P <- df$DDC__Cardiometabolic_o2csf__NPX
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_2 <- lines.roc(obs, pred,  col="#E69F00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)

N_P <- df$CSF_NFL_pgml_Imputed_NTK_2020
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_3 <- lines.roc(obs, pred,  col="#009E73", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)

N_P <- df$`CRH-DDC`
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_4 <- lines.roc(obs, pred,  col="#CC79A7", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)

N_P <- df$`CRH-NFL`
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_5 <- lines.roc(obs, pred,  col="#0072B2", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)

N_P <- df$`CRH+NFL-DDC`
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_6 <- lines.roc(obs, pred,  col="#56B4E9", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)

N_P <- df$`CRH-DDC+NFL`
model <- glm (Diagnosis ~ N_P + age + gender_baseline_variable, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
obs = model$y
par(cex.axis=1.5)
p_7 <- lines.roc(obs, pred,  col="#D55E00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2)
