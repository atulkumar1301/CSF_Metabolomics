#Box Plot
library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("Full_Data_Metabolites.txt")
df$aSyn_SAA_Status <- with(df, ifelse(df$SAA_Status == 1, "aSyn-SAA+", "aSyn-SAA-"))
av <- aov(`carboxyethyl-GABA` ~ aSyn_SAA_Status, data = df)
summary (av)
my_comparisons <- list( c("aSyn-SAA+", "aSyn-SAA-"))

p <- ggplot(data=subset(df, !is.na(aSyn_SAA_Status)), aes(x = aSyn_SAA_Status, y = `3-methoxytyrosine`,
                    colour = aSyn_SAA_Status,
                    shape = aSyn_SAA_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p <- p + xlab (expression(alpha*"Syn SAA Status")) + ylab (expression ("3-Methoxytyrosine"))
p <- p + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("a) 3-Methoxytyrosine"))
p

p_1 <- ggplot(data=subset(df, !is.na(aSyn_SAA_Status)), aes(x = aSyn_SAA_Status, y = `dopamine 3-O-sulfate`,
                      colour = aSyn_SAA_Status,
                      shape = aSyn_SAA_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_1 <- p_1 + xlab (expression(alpha*"Syn SAA Status")) + ylab (expression ("Dopamine 3-O-sulfate"))
p_1 <- p_1 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("b) Dopamine 3-O-sulfate"))
p_1

p_2 <- ggplot(data=subset(df, !is.na(aSyn_SAA_Status)), aes(x = aSyn_SAA_Status, y = `3-methoxytyramine sulfate`,
                      colour = aSyn_SAA_Status,
                      shape = aSyn_SAA_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_2 <- p_2 + xlab (expression(alpha*"Syn SAA Status")) + ylab (expression ("3-methoxytyramine sulfate"))
p_2 <- p_2 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("c) 3-methoxytyramine sulfate"))
p_2

p_3 <- ggplot(data=subset(df, !is.na(aSyn_SAA_Status)), aes(x = aSyn_SAA_Status, y = `3-hydroxyhexanoate`,
                      colour = aSyn_SAA_Status,
                      shape = aSyn_SAA_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_3 <- p_3 + xlab (expression(alpha*"Syn SAA Status")) + ylab (expression ("3-Hydroxyhexanoate"))
p_3 <- p_3 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("d) 3-Hydroxyhexanoate"))
p_3

p_4 <- ggplot(data=subset(df, !is.na(aSyn_SAA_Status)), aes(x = aSyn_SAA_Status, y = `N-acetylhistidine`,
                      colour = aSyn_SAA_Status,
                      shape = aSyn_SAA_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_4 <- p_4 + xlab (expression(alpha*"Syn SAA Status")) + ylab (expression ("N-Acetylhistidine"))
p_4 <- p_4 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("e) N-Acetylhistidine"))
p_4


plot <- ggarrange(p, p_1, p_2, p_3, p_4, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Top five metabolites associated with "*alpha* "Syn SAA Status"), color = "#999999", face = "bold", size = 20, family = "serif"))
