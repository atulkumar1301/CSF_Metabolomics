#Box Plot
library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("Full_Data_Metabolites.txt")
df$WML_Status <- with(df, ifelse(df$samseg_icv_tertiles == 1, "WML+", "WML-"))
av <- aov(`carboxyethyl-GABA` ~ WML_Status, data = df)
summary (av)
my_comparisons <- list( c("WML+", "WML-"))

p <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = `1-(1-enyl-stearoyl)-2-docosahexaenoyl-GPE (P-18:0/22:6)*`,
                    colour = WML_Status,
                    shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p <- p + xlab (expression("WML Status")) + ylab (expression ("1-(1-Enyl-stearoyl)-2-docosahexaenoyl-gpe"))
p <- p + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("a) 1-(1-Enyl-stearoyl)-2-docosahexaenoyl-gpe"))
p

p_1 <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = `N-acetyl-aspartyl-glutamate (NAAG)`,
                      colour = WML_Status,
                      shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_1 <- p_1 + xlab (expression("WML Status")) + ylab (expression ("N-acetyl-aspartyl-glutamate"))
p_1 <- p_1 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("b) N-acetyl-aspartyl-glutamate"))
p_1

p_2 <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = `X-10457`,
                      colour = WML_Status,
                      shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_2 <- p_2 + xlab (expression("WML Status")) + ylab (expression ("X-10457"))
p_2 <- p_2 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("c) X-10457"))
p_2

p_3 <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = `5-methylthioribose**`,
                      colour = WML_Status,
                      shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_3 <- p_3 + xlab (expression("WML Status")) + ylab (expression ("5−Methylthioribose"))
p_3 <- p_3 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("d) 5−Methylthioribose"))
p_3

p_4 <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = allantoin,
                      colour = WML_Status,
                      shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_4 <- p_4 + xlab (expression("WML Status")) + ylab (expression ("Allantoin"))
p_4 <- p_4 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("e) Allantoin"))
p_4

p_5 <- ggplot(data=subset(df, !is.na(WML_Status)), aes(x = WML_Status, y = `phenyllactate (PLA)`,
                      colour = WML_Status,
                      shape = WML_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_5 <- p_5 + xlab (expression("WML Status")) + ylab (expression ("Phenyllactate"))
p_5 <- p_5 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=12),
        legend.text = element_text(family = "serif", size=12),
        panel.background = element_blank()) + labs(title=expression("f) Phenyllactate"))
p_5


plot <- ggarrange(p, p_1, p_2, p_3, p_4, p_5, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Top six metabolites associated with WML Status"), color = "#999999", face = "bold", size = 20, family = "serif"))
