#Box Plot
library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("Full_Data_Metabolites.txt")
df$Ab_Status <- with(df, ifelse(df$Abnormal_CSF_Ab42_Ab40_Ratio == 1, "Ab+", "Ab-"))
av <- aov(`carboxyethyl-GABA` ~ Ab_Status, data = df)
summary (av)
my_comparisons <- list( c("Ab+", "Ab-"))

p <- ggplot(df, aes(x = Ab_Status, y = `X-25790`,
                    colour = Ab_Status,
                    shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p <- p + xlab (expression("A"*beta*" Status")) + ylab (expression ("X-25790"))
p <- p + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("a) X-25790"))
p

p_1 <- ggplot(df, aes(x = Ab_Status, y = df$`X-18887`,
                    colour = Ab_Status,
                    shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_1 <- p_1 + xlab (expression("A"*beta*" Status")) + ylab (expression ("X-18887"))
p_1 <- p_1 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("b) X-18887"))
p_1

p_2 <- ggplot(df, aes(x = Ab_Status, y = `carboxyethyl-GABA`,
               colour = Ab_Status,
               shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_2 <- p_2 + xlab (expression("A"*beta*" Status")) + ylab (expression ("Carboxyethyl-GABA"))
p_2 <- p_2 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("c) Carboxyethyl-GABA"))
p_2

p_3 <- ggplot(df, aes(x = Ab_Status, y = `3-methyl-2-oxobutyrate`,
                      colour = Ab_Status,
                      shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_3 <- p_3 + xlab (expression("A"*beta*" Status")) + ylab (expression ("3-Methyl-2-oxobutyrate"))
p_3 <- p_3 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("d) 3-Methyl-2-oxobutyrate"))
p_3

p_4 <- ggplot(df, aes(x = Ab_Status, y = `gamma-glutamylleucine`,
                      colour = Ab_Status,
                      shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_4 <- p_4 + xlab (expression("A"*beta*" Status")) + ylab (expression ("gamma−Glutamyl−leucine"))
p_4 <- p_4 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("e) gamma−Glutamyl−leucine"))
p_4

p_5 <- ggplot(df, aes(x = Ab_Status, y = `4-methyl-2-oxopentanoate`,
                      colour = Ab_Status,
                      shape = Ab_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p_5 <- p_5 + xlab (expression("A"*beta*" Status")) + ylab (expression ("4-Methyl-2-oxopentanoate"))
p_5 <- p_5 + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("f) 4-Methyl-2-oxopentanoate"))
p_5


plot <- ggarrange(p, p_1, p_2, p_3, p_4, p_5, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Top six metabolites associated with A"*beta*" Status"), color = "#999999", face = "bold", size = 20, family = "serif"))
