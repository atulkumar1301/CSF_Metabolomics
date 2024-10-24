#Plotting ROC Curve
library(data.table)
library(ROCR)
library(pROC)

df <- fread ("Full_Data_Metabolites.txt")
par(mfrow=c(2,2))

###Abeta
model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `X-25790` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot(roc(obs, pred),  col="#F0E442", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.8, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, main = expression("a) ROC Curve: Metabolites as a predictor for A"*beta* " Status"), cex.main = 1)
legend("bottomright", legend=c("X-25790", "X-18887", "Carboxyethyl−GABA", "3−Methyl−2−oxobutyrate", "gamma−Glutamyl−leucine", "4−Methyl−2−oxopentanoate"), col=c("#F0E442", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"), lwd=2, cex = 0.75)


model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `X-18887` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot (roc(obs, pred),  col="#E69F00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.75, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `carboxyethyl-GABA` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

par(cex.axis=1)
plot (roc(obs, pred),  col="#009E73", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.7, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `3-methyl-2-oxobutyrate` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#CC79A7", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.65, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `gamma-glutamylleucine` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#0072B2", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.6, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (Abnormal_CSF_Ab42_Ab40_Ratio ~ `4-methyl-2-oxopentanoate` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#D55E00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.6, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)

####Tau

model <- glm (Tau_Binary ~ `carboxyethyl-GABA` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + Abnormal_CSF_Ab42_Ab40_Ratio + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot(roc(obs, pred),  col="#F0E442", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.8, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, main = expression("b) ROC Curve: Metabolites as a predictor for TauPET Status"), cex.main = 1)
legend("bottomright", legend=c("Carboxyethyl−GABA"), col=c("#F0E442"), lwd=2, cex = 0.75)


###SAA

model <- glm (SAA_Status ~ `3-methoxytyrosine` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot(roc(obs, pred),  col="#F0E442", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.8, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, main = expression("c) ROC Curve: Metabolites as a predictor for "*alpha*"Syn SAA Status"), cex.main = 1)
legend("bottomright", legend=c("3−Methoxytyrosine", "Dopamine 3−O−sulfate", "3−methoxytyramine sulfate", "3−Hydroxyhexanoate", "N−Acetylhistidine"), col=c("#F0E442", "#E69F00", "#009E73", "#56B4E9", "#D55E00"), lwd=2, cex = 0.75)


model <- glm (SAA_Status ~ `dopamine 3-O-sulfate` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot (roc(obs, pred),  col="#E69F00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.75, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (SAA_Status ~ `3-methoxytyramine sulfate` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

par(cex.axis=1)
plot (roc(obs, pred),  col="#009E73", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.7, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (SAA_Status ~ `3-hydroxyhexanoate` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#0072B2", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.65, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (SAA_Status ~ `N-acetylhistidine` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + Abnormal_CSF_Ab42_Ab40_Ratio + tnic_cho_com_I_IV + samseg_wmhs_WMH_total_mm3, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#D55E00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.6, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


###WML
model <- glm (samseg_icv_tertiles ~ `1-(1-enyl-palmitoyl)-2-oleoyl-GPC (P-16:0/18:1)*` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot(roc(obs, pred),  col="#F0E442", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.8, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, main = expression("d) ROC Curve: Metabolites as a predictor for WML Status"), cex.main = 1)
legend("bottomright", legend=c("1−(1−Enyl−stearoyl)−2−docosahexaenoyl−gpe", "N−acetyl−aspartyl−glutamate", "X−10457", "5−Methylthioribose", "Allantoin", "Phenyllactate"), col=c("#F0E442", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"), lwd=2, cex = 0.75)


model <- glm (samseg_icv_tertiles ~ `N-acetyl-aspartyl-glutamate (NAAG)` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")


par(cex.axis=1)
plot (roc(obs, pred),  col="#E69F00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.75, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (samseg_icv_tertiles ~ `X-10457` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

par(cex.axis=1)
plot (roc(obs, pred),  col="#009E73", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.7, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (samseg_icv_tertiles ~ `5-methylthioribose**` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#CC79A7", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.65, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (samseg_icv_tertiles ~ allantoin + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#0072B2", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.6, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)


model <- glm (samseg_icv_tertiles ~ `phenyllactate (PLA)` + Age + Gender + mean_standardized_metabolomic_level + Recruitment_Bias + SAA_Status + tnic_cho_com_I_IV + Abnormal_CSF_Ab42_Ab40_Ratio, data = df, family = binomial (link = 'logit'))
obs = model$y
pred = model$fitted.values
auc (obs, pred)
ci.auc (obs, pred)
aa <- roc (obs, pred)
best.coords <- coords(aa, "best", best.method="youden")

plot (roc(obs, pred),  col="#D55E00", col.axis="black", col.lab="black", lwd=2, asp = NA, indentity = TRUE, identity.lty = "dashed", identity.col = "black", identity.lwd = 2, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.6, print.auc.pattern = "AUC = %.2f", print.auc.cex=1, cex.lab = 1, cex.main = 1, add = TRUE)

