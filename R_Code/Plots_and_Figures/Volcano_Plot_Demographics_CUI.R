#Volcano plot
library(ggpubr)
library(data.table)
library(ggplot2)
library(ggrepel)
library (plotly)


# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#D55E00", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#F0E442", "#CC79A7")

####Age

df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/7_Full_Result_Data_Analysis_Age.txt")
p <- ggplot (data = df, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.01768617021), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00518617021 (with No Dynamics)
p <- p + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p <- p + scale_x_continuous(breaks = round(seq(-0.08, 0.06, by = 0.01),1))
p <- p + scale_y_continuous(breaks = round (seq (0, 18, by = 2), 1))
p <- p + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p <- p +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("a) Age"))
p


###Gender

df_2 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/TauPET/7_Full_Result_Data_Analysis_Taupet.txt")
p_2 <- ggplot (data = df_2, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_2 <- p_2 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.02087765957), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00119680851 (with No Dynamics)
p_2 <- p_2 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_2 <- p_2 + scale_x_continuous(breaks = round(seq(-0.7, 0.7, by = 0.2),1))
p_2 <- p_2 + scale_y_continuous(breaks = round (seq (0, 23, by = 2), 1))
p_2 <- p_2 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_2 <- p_2 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("b) Sex"))
p_2



plot <- ggarrange(p, p_2, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Individual effects of age and sex on CSF Metabolites in CUI Individuals"), color = "#999999", face = "bold", size = 20, family = "serif"))
