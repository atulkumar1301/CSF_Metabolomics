####Age Metabolites

df <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/Full/7_Full_Result_Data_Analysis_Age_Full.txt")
p <- ggplot (data = df, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.02087765957), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00518617021 (with No Dynamics)
p <- p + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p <- p + scale_x_continuous(breaks = round(seq(-0.07, 0.04, by = 0.01),1))
p <- p + scale_y_continuous(breaks = round (seq (0, 30, by = 2), 1))
p <- p + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p <- p +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("A Age: Metabolites"))
p


###Gender Metabolites

df_2 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/Full/7_Full_Result_Data_Analysis_Gender_Full.txt")
p_2 <- ggplot (data = df_2, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_2 <- p_2 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.03271276595), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00119680851 (with No Dynamics)
p_2 <- p_2 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_2 <- p_2 + scale_x_continuous(breaks = round(seq(-0.8, 0.6, by = 0.2),1))
p_2 <- p_2 + scale_y_continuous(breaks = round (seq (0, 55, by = 4), 1))
p_2 <- p_2 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_2 <- p_2 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("B Sex: Metabolites"))
p_2

####CU_CI Metabolites

df_4 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Demographics/Full/7_Full_Result_Data_Analysis_CU_CI_Full.txt")
p_4 <- ggplot (data = df_4, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_4 <- p_4 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.01276595744), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00864361702 (with No Dynamics)
p_4 <- p_4 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_4 <- p_4 + scale_x_continuous(breaks = round(seq(-0.8, 0.5, by = 0.2), 1))
p_4 <- p_4 + scale_y_continuous(breaks = round (seq (0, 15, by = 2), 1))
p_4 <- p_4 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_4 <- p_4 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("C CU/CI Status Metabolites"))
p_4


####Age Sub-Pathways

df_6 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/Full/7_Full_Result_Data_Analysis_Subpathway_Age_Full.txt")
p_6 <- ggplot (data = df_6, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_6 <- p_6 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p_6-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.02183098591), linetype = "FDR p_6-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p_6-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00518617021 (with No Dynamics)
p_6 <- p_6 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_6 <- p_6 + scale_x_continuous(breaks = round(seq(-0.02, 0.02, by = 0.01),1))
p_6 <- p_6 + scale_y_continuous(breaks = round (seq (0, 30, by = 2), 1))
p_6 <- p_6 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(P)))
p_6 <- p_6 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("D Age: Sub-Pathways"))
p_6


###Gender Sub-Pathways

df_7 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/Full/7_Full_Result_Data_Analysis_Subpathway_Gender_Full.txt")
p_7 <- ggplot (data = df_7, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_7 <- p_7 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.03591549295), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00119680851 (with No Dynamics)
p_7 <- p_7 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_7 <- p_7 + scale_x_continuous(breaks = round(seq(-0.5, 0.4, by = 0.2),1))
p_7 <- p_7 + scale_y_continuous(breaks = round (seq (0, 55, by = 3), 1))
p_7 <- p_7 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_7 <- p_7 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("E Sex: Sub-Pathways"))
p_7

####CU_CI Sub-Pathways

df_8 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Demographics/Full/7_Full_Result_Data_Analysis_Subpathway_CU_CI_Full.txt")
p_8 <- ggplot (data = df_8, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p_8 <- p_8 + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.01056338028), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00")))) ##0.00864361702 (with No Dynamics)
p_8 <- p_8 + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p_8 <- p_8 + scale_x_continuous(breaks = round(seq(-0.2, 0.3, by = 0.2), 1))
p_8 <- p_8 + scale_y_continuous(breaks = round (seq (0, 15, by = 2), 1))
p_8 <- p_8 + xlab ("Effect") + labs (color = "Regulation") + ylab (expression (-log[10]~(p)))
p_8 <- p_8 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=12),
        axis.text.y = element_text(family = "serif", size=12),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title=expression("F CU/CI Status Sub-Pathways"))
p_8


plot <- ggarrange(p, p_2, p_4, p_6, p_7, p_8, common.legend = TRUE, legend = "bottom")

plot
#annotate_figure(plot, top = text_grob(expression("Effects of demographics on CSF Metabolites"), color = "#999999", face = "bold", size = 20, family = "serif"))
