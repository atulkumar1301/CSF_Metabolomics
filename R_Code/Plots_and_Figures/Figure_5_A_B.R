library (ggplot2)
library (reshape2)
library(data.table)
library(ggpubr)

#Not Adjusted
df_1 <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Heat_Map_Data.txt", header = TRUE, stringsAsFactors = FALSE)
ord <- hclust( dist(df_1, method = "euclidean"), method = "ward.D" )$order
melted_df <- melt (df)
ggheatmap <- ggplot (data = melted_df, aes (x = Metabolite, y = variable, fill = value)) + 
  geom_tile (color = "white") +
  scale_x_discrete(limits=melted_df$Metabolite[ord]) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.044276440, limit = c(-5, 6), space = "Lab", name="Effect Size") +
  xlab("Metabolites") +
  ylab("Pathologies") +
  ggtitle ("A")

ggheatmap <- ggheatmap +
  theme(
    plot.title = element_text(family = "serif", size=14, face = "bold"),
    axis.title.x = element_text(family = "serif", size=12),
    axis.title.y = element_text(family = "serif", size=12),
    axis.text.x = element_text(family = "serif", size=9, angle = 90),
    axis.text.y = element_text(family = "serif", size=12),
    legend.title = element_text(family = "serif", size=12),
    legend.text = element_text(family = "serif", size=12),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggheatmap


#Adjusted
df <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Heat_Map_Data_Adjusted.txt")
ord <- hclust( dist(df_1, method = "euclidean"), method = "ward.D" )$order
melted_df <- melt (df)
ggheatmap_1 <- ggplot (data = melted_df, aes (x = Metabolite, y = variable, fill = value)) + 
  geom_tile (color = "white") +
  scale_x_discrete(limits=melted_df$Metabolite[ord]) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.110758840, limit = c(-5, 6), space = "Lab", name="Effect Size")+
  xlab("Metabolites") +
  ylab("Pathologies") +
  ggtitle ("B Adjusted for other co-pathologies")
ggheatmap_1 <- ggheatmap_1 +
  theme(
    plot.title = element_text(family = "serif", size=14, face = "bold"),
    axis.title.x = element_text(family = "serif", size=12),
    axis.title.y = element_text(family = "serif", size=12),
    axis.text.x = element_text(family = "serif", size=9, angle = 90),
    axis.text.y = element_text(family = "serif", size=12),
    legend.title = element_text(family = "serif", size=12),
    legend.text = element_text(family = "serif", size=12),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggheatmap_1

plot <- ggarrange(ggheatmap, ggheatmap_1, common.legend = TRUE, legend = "left", ncol = 1, nrow = 2)

plot
