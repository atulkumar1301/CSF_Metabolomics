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


####Age adjusted for Gender and CU/CI: Metabolites

df_1 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/All_Combined/7_Full_Result_Data_Analysis_Age_All.txt")
p_1 <- ggplot (data = df_1, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_1 <- p_1 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.02180851063), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00039893617 (with no Dynamics)
p_1 <- p_1 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_1 <- p_1 + scale_x_continuous(breaks = round(seq(-0.08, 0.07, by = 0.01),1))
p_1 <- p_1 + scale_y_continuous(breaks = round (seq (0, 24, by = 2), 1))
p_1 <- p_1 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p_1 <- p_1 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("A Age adjusted for Sex and CU/CI status: Metabolites"))
p_1


###Gender adjusted for Age and CU/CI: Metabolites

df_3 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/All_Combined/7_Full_Result_Data_Analysis_Gender_All.txt")
p_3 <- ggplot (data = df_3, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_3 <- p_3 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.0222074468), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00053191489 (with No Dynamics)
p_3 <- p_3 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_3 <- p_3 + scale_x_continuous(breaks = round(seq(-0.7, 0.6, by = 0.2),1))
p_3 <- p_3 + scale_y_continuous(breaks = round (seq (0, 58, by = 4), 1))
p_3 <- p_3 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_3 <- p_3 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("B Sex adjusted for Age and CU/CI status: Metabolites"))
p_3


##CU/CI adjusted for Gender and Age: Metabolites

df_5 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/All_Combined/7_Full_Result_Data_Analysis_CU_CI_All.txt")
p_5 <- ggplot (data = df_5, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_5 <- p_5 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00851063829), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.01316489361 (with No Dynamics)
p_5 <- p_5 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_5 <- p_5 + scale_x_continuous(breaks = round(seq(-0.3, 0.4, by = 0.2), 1))
p_5 <- p_5 + scale_y_continuous(breaks = round (seq (0, 17, by = 1), 1))
p_5 <- p_5 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_5 <- p_5 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("C CU/CI status adjusted for Age and Sex: Metabolites"))
p_5



####Age adjusted for Gender and CU/CI Sub-Pathways

df_6 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/All_Combined/7_Full_Result_Data_Analysis_Subpathway_Age_All_Combined.txt")
p_6 <- ggplot (data = df_6, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_6 <- p_6 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.02042253521), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00039893617 (with no Dynamics)
p_6 <- p_6 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_6 <- p_6 + scale_x_continuous(breaks = round(seq(-0.02, 0.02, by = 0.01),1))
p_6 <- p_6 + scale_y_continuous(breaks = round (seq (0, 24, by = 1), 1))
p_6 <- p_6 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p_6 <- p_6 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("D Age adjusted for Sex and CU/CI status: Sub-Pathways"))
p_6


###Gender adjusted for Age and CU/CI Sub-Pathways

df_7 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/All_Combined/7_Full_Result_Data_Analysis_Subpathway_Gender_All_Combined.txt")
p_7 <- ggplot (data = df_7, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_7 <- p_7 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.0222074468), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00053191489 (with No Dynamics)
p_7 <- p_7 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_7 <- p_7 + scale_x_continuous(breaks = round(seq(-0.5, 0.4, by = 0.2),1))
p_7 <- p_7 + scale_y_continuous(breaks = round (seq (0, 58, by = 3), 1))
p_7 <- p_7 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_7 <- p_7 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("E Sex adjusted for Age and CU/CI status: Sub-Pathways"))
p_7


##CU/CI adjusted for Gender and Age Sub-Pathways

df_8 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/All_Combined/7_Full_Result_Data_Analysis_Subpathway_CU_CI_All_Combined.txt")
p_8 <- ggplot (data = df_8, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_8 <- p_8 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00704225352), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.01316489361 (with No Dynamics)
p_8 <- p_8 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_8 <- p_8 + scale_x_continuous(breaks = round(seq(-0.2, 0.3, by = 0.1), 1))
p_8 <- p_8 + scale_y_continuous(breaks = round (seq (0, 17, by = 1), 1))
p_8 <- p_8 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_8 <- p_8 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=11, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("F CU/CI status adjusted for Age and Sex: Sub-Pathways"))
p_8



plot <- ggarrange(p_1, p_3, p_5, p_6, p_7, p_8, common.legend = TRUE, legend = "bottom")

plot
#annotate_figure(plot, top = text_grob(expression("Effects of demographics on CSF Metabolites"), color = "#999999", face = "bold", size = 20, family = "serif"))




