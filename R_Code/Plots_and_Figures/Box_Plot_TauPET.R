#Box Plot
library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("Full_Data_Metabolites.txt")
df$TauPET_Status <- with(df, ifelse(df$Tau_Binary == 1, "TauPET+", "TauPET-"))
av <- aov(`carboxyethyl-GABA` ~ TauPET_Status, data = df)
summary (av)
my_comparisons <- list( c("TauPET+", "TauPET-"))

p <- ggplot(data=subset(df, !is.na(TauPET_Status)), aes(x = TauPET_Status, y = `carboxyethyl-GABA`,
                    colour = TauPET_Status,
                    shape = TauPET_Status)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  geom_jitter() + theme_gray() +
  stat_compare_means(comparisons = my_comparisons)
p <- p + xlab (expression("TauPET Status")) + ylab (expression ("Carboxyethyl-GABA"))
p <- p + theme_light()+ 
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("Carboxyethyl-GABA"))
p


plot <- ggarrange(p, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Top metabolite associated with TauPET Status"), color = "#999999", face = "bold", size = 20, family = "serif"))
