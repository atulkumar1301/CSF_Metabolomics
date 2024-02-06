#Volcano plot
library(ggpubr)
library(data.table)
library(ggplot2)
library(ggrepel)
library (plotly)
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/Age_Stratified/Age_Stratified_Greater_Than_60/4_Full_Result_Data_Analysis_Raw_Adjusted.txt")
p <- ggplot (data = df, aes (x = Effect, y = -log10(P), col = Regulation, label = Label))+
  geom_point ()+ geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_hline (aes(yintercept=-log10(0.05), linetype = "p-value 0.05", col="black")) +
  geom_hline (aes (yintercept=-log10(0.014), linetype = "FDR p-value 0.05", col="#D55E00")) +
  scale_linetype_manual(name = "p-value cut off", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("black", "#D55E00"))))
#p <- p + geom_vline(xintercept=c(-0.165879033,  0.128549555), col="black")
p <- p + scale_color_manual(values=cbbPalette, limits = force) + theme_light()
p <- p + scale_x_continuous(breaks = round(seq(-2.9, 5.6, by = 1),1))
p <- p + scale_y_continuous(breaks = round (seq (0, 12, by = 1), 1))
p <- p + xlab ("Effect") + labs (color = "Regulation")
p <- p +
  theme(legend.position="left",
               plot.title = element_text(family = "serif", size=18, face = "bold", hjust = 0.5),
               axis.title.x = element_text(family = "serif", size=16),
               axis.title.y = element_text(family = "serif", size=16),
               axis.text.x = element_text(family = "serif", size=10, angle = 10),
               axis.text.y = element_text(family = "serif", size=10),
               legend.title = element_text(family = "serif", size=16),
               legend.text = element_text(family = "serif", size=16),
               panel.background = element_blank()) + labs(title="Differential Regulation of CSF Metabolites with respect to Sex \n for Age > 60 Years")# + labs(title=expression("Atypical vs A"*beta*"+ MCI and AD"))
p
