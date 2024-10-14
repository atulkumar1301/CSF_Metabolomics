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

####Abeta

df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/A_Beta/Medication/7_Full_Result_Data_Analysis_A_Beta_Medication.txt")
p <- ggplot (data = df, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00664893617), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00518617021 (with No Dynamics)
p <- p + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p <- p + scale_x_continuous(breaks = round(seq(-0.9, 0.3, by = 0.3),1))
p <- p + scale_y_continuous(breaks = round (seq (0, 10, by = 1), 1))
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
        panel.background = element_blank()) + labs(title=expression("a) A"*beta))
p


###TauPET

df_2 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/TauPET/Medication/7_Full_Result_Data_Analysis_Taupet_Medication.txt")
p_2 <- ggplot (data = df_2, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_2 <- p_2 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00106382978), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00119680851 (with No Dynamics)
p_2 <- p_2 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_2 <- p_2 + scale_x_continuous(breaks = round(seq(-0.5, 0.4, by = 0.2),1))
p_2 <- p_2 + scale_y_continuous(breaks = round (seq (0, 11, by = 1), 1))
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
        panel.background = element_blank()) + labs(title=expression("b) TauPET"))
p_2

####Asyn

df_4 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/A_Syn/Medication/7_Full_Result_Data_Analysis_A_Syn_Medication.txt")
p_4 <- ggplot (data = df_4, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_4 <- p_4 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00186170212), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00864361702 (with No Dynamics)
p_4 <- p_4 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_4 <- p_4 + scale_x_continuous(breaks = round(seq(-0.5, 0.8, by = 0.2), 1))
p_4 <- p_4 + scale_y_continuous(breaks = round (seq (0, 15, by = 2), 1))
p_4 <- p_4 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_4 <- p_4 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("c) SAA"))
p_4

####WML

df_6 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/WML/Medication/7_Full_Result_Data_Analysis_WML_Medication.txt")
p_6 <- ggplot (data = df_6, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_6 <- p_6 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00864361702), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00784574468 (with No Dynamics)
p_6 <- p_6 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_6 <- p_6 + scale_x_continuous(breaks = round(seq(-0.9, 1, by = 0.3), 1))
p_6 <- p_6 + scale_y_continuous(breaks = round (seq (0, 8, by = 1), 1))
p_6 <- p_6 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_6 <- p_6 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("d) WML"))
p_6

plot <- ggarrange(p, p_2, p_4, p_6, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Individual effects of all pathologies on CSF Metabolites adjusted for Medication"), color = "#999999", face = "bold", size = 20, family = "serif"))


####Abeta adjusted for other pathologies

df_1 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/7_Full_Result_Data_Analysis_Ab.txt")
p_1 <- ggplot (data = df_1, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_1 <- p_1 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00319148936), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00039893617 (with no Dynamics)
p_1 <- p_1 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_1 <- p_1 + scale_x_continuous(breaks = round(seq(-0.3, 0.3, by = 0.1),1))
p_1 <- p_1 + scale_y_continuous(breaks = round (seq (0, 8, by = 1), 1))
p_1 <- p_1 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p_1 <- p_1 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("a) A"*beta))
p_1


###TauPET adjusted for other pathologies

df_3 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/7_Full_Result_Data_Analysis_Taupet.txt")
p_3 <- ggplot (data = df_3, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_3 <- p_3 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00026595744), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00053191489 (with No Dynamics)
p_3 <- p_3 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_3 <- p_3 + scale_x_continuous(breaks = round(seq(-0.6, 0.4, by = 0.2),1))
p_3 <- p_3 + scale_y_continuous(breaks = round (seq (0, 8, by = 1), 1))
p_3 <- p_3 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_3 <- p_3 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("b) TauPET"))
p_3


##Asyn adjusted for other pathologies

df_5 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/7_Full_Result_Data_Analysis_Asyn.txt")
p_5 <- ggplot (data = df_5, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_5 <- p_5 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00146276595), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.01316489361 (with No Dynamics)
p_5 <- p_5 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_5 <- p_5 + scale_x_continuous(breaks = round(seq(-0.3, 0.7, by = 0.3), 1))
p_5 <- p_5 + scale_y_continuous(breaks = round (seq (0, 17, by = 1), 1))
p_5 <- p_5 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_5 <- p_5 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("c)SAA"))
p_5


###WML adjusted for other pathologies

df_7 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/All_Combined/7_Full_Result_Data_Analysis_WML.txt")
p_7 <- ggplot (data = df_7, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_7 <- p_7 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.00917553191), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.0097074468 (with No Dynamics)
p_7 <- p_7 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_7 <- p_7 + scale_x_continuous(breaks = round(seq(-1, 1, by = 0.3), 1))
p_7 <- p_7 + scale_y_continuous(breaks = round (seq (0, 15, by = 1), 1))
p_7 <- p_7 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_7 <- p_7 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("d) WML"))
p_7


plot <- ggarrange(p_1, p_3, p_5, p_7, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob(expression("Effects of all neurodegenerative disease pathology on CSF Metabolites when adjusted for other pathologies"), color = "#999999", face = "bold", size = 20, family = "serif"))
