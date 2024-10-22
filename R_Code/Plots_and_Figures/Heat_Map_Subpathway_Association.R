library (ggplot2)
library (reshape2)
library(data.table)
library(ggpubr)
#Including APOE Region
df_1 <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Heat_Map_Data.txt", header = TRUE, stringsAsFactors = FALSE)
ord <- hclust( dist(df_1, method = "euclidean"), method = "ward.D" )$order
melted_df <- melt (df)
ggheatmap <- ggplot (data = melted_df, aes (x = Metabolite, y = variable, fill = value)) + 
  geom_tile (color = "white") +
  scale_x_discrete(limits=melted_df$Metabolite[ord]) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.044276440, limit = c(-5, 6), space = "Lab", name="Effect Size") +
  xlab("Metabolites") +
  ylab("Pathologies") +
  ggtitle ("Metabolite association with Different Pathologies")

ggheatmap <- ggheatmap +
  theme(
    plot.title = element_text(family = "serif", size=18, face = "bold"),
    axis.title.x = element_text(family = "serif", size=16),
    axis.title.y = element_text(family = "serif", size=16),
    axis.text.x = element_text(family = "serif", size=9, angle = 90),
    axis.text.y = element_text(family = "serif", size=14),
    legend.title = element_text(family = "serif", size=16),
    legend.text = element_text(family = "serif", size=16),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggheatmap
