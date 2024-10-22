#Box Plot
library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("Full_Data_Metabolites.txt")
df$Ab <- with(df, ifelse(df$Abnormal_CSF_Ab42_Ab40_Ratio == 1, "Ab+", "Ab-"))
av <- aov(`carboxyethyl-GABA` ~ Ab, data = df)
summary (av)
my_comparisons <- list( c("Ab+", "Ab-"))
ggplot(df, aes(x = Ab, y = `carboxyethyl-GABA`,
               colour = Ab,
               shape = Ab)) + 
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
stat_compare_means(comparisons = my_comparisons)
#p <- p + scale_y_continuous(breaks=seq(-5, 0.5, 2.5), limits = c(-5, 2.5))
p <- p + stat_compare_means(aes (label = paste0("p = ", after_stat(p.format))), label.x = 2.5, label.y = -5)
#p <- p + stat_compare_means(aes (label = paste0("p", after_stat(p.format))))
p <- p + stat_compare_means(comparisons = my_comparisons)
#p <- p + theme(axis.text.x = element_text(angle = 15))
#p <- p + theme(text = element_text(family = "Times New Roman"))
dev.off()
p
ggplot(df, aes(x = Group, y = DDC,
               colour = Group,
               shape = Group)) + 
  geom_jitter()
my_comparisons <- list( c("CTR", "De Novo LBD"))
p_1 <- ggplot(df, aes(x = Group, y = DDC, colour = Group)) + geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_color_manual(values = c("black", "red")) +
  geom_jitter() + theme_gray() +
  theme(legend.position="left") +
  labs(title="B) CTR vs De Novo LBD", face = "bold") +
  theme(legend.position="none") +
  ylab ("Plasma DDC") +
  theme(axis.title.x=element_blank()) +
  stat_compare_means(comparisons = my_comparisons)
p_1



plot <- ggarrange(p_1, p_2, p_3, p_4, p_5, p_6, ncol = 3, nrow = 2)

annotate_figure(plot, top = text_grob("BioFINDER-1 cohort (Plasma)", color = "red", face = "bold", size = 18))

library(data.table)
par(mfrow=c(1,2)) # for multiple plots
#par (bg = "lightgrey")
df <- fread ("Full_Data_Plasma_BF-1.txt")
par(cex.axis = 1.5)
boxplot(DDC ~ Group, data = df, lwd = 2, ylab = 'Plasma DDC', notch = TRUE, cex.lab = 1.5, main = "A) Box plot: CTR vs LBD vs Atypical PS", cex.main = 1.5)
stripchart(DDC ~ Group, vertical = TRUE, data = df, 
           method = "jitter", add = TRUE, pch = 20, col = c('black', 'red', 'green'), cex = 2, jitter = 0.3)

mtext("BioFINDER-1 Cohort (Plasma)", side = 3, line = -2, outer = TRUE, cex = 2, col = 'red')
df_1 <- fread ("Data_SAA-_CTR_vs_SAA+_LBD.txt")
par(cex.axis=1.5)
boxplot(DDC ~ Group, data = df_1, lwd = 2, ylab = 'CSF DDC', notch = TRUE, cex.lab = 1.5, main = "B) SAA- CTR vs SAA+ LBD")
stripchart(DDC ~ Group, vertical = TRUE, data = df, 
           method = "jitter", add = TRUE, pch = 20, col = c('black', 'red'))

df_2 <- fread ("Data_SAA-_CTR_vs_SAA+_LBD.txt")
par(cex.axis=1.5)
df_2
boxplot(DDC ~ Group, data = df_1, lwd = 2, ylab = 'CSF DDC', notch = TRUE, cex.lab = 1.5, main = "B) SAA- CTR vs SAA+ LBD")
stripchart(DDC ~ Group, vertical = TRUE, data = df, 
           method = "jitter", add = TRUE, pch = 20, col = c('black', 'red'))






