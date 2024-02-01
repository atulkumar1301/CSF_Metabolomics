library(data.table)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(factoextra)
library(ggplot2)
library(ggrepel)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/5_Data_Analysis_1_Full_Imputed.txt")

### PCA
df_PCA <- prcomp (df [,c (30:398)], center = TRUE, scale. = TRUE)
fviz_eig(df_PCA)

####Extract Loading Data
df_Rotation <- as.data.frame(df_PCA$rotation)
write.table(df_Rotation, file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/9_Data_Analysis_1_PCA_Rotation.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

###Loading Plot
metabolite_selected <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/10_Data_Analysis_1_PCA_Loading.txt")

p <- ggplot(metabolite_selected, aes(x = PC1, y = PC2, label = Metabolites), size=3, vjust="outward")+
  geom_point() + geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, color="red")
p <- p + theme_bw()
p <- p +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold"),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title="Loadings Plot")
p
