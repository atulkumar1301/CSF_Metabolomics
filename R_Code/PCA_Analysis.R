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

######Find the number for PC
###The variance is the square of the standard deviation
variance <- df_PCA$sdev^2

###The proportion of variance is the variance divided by the sum of all variances
prop_variance <- variance / sum (variance)

###Cumulative Proportion is the cumulative sum of the proportion of variance
cum_Prop <- cumsum (prop_variance)

write.table(cum_Prop, file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/12_Cumulative_Variance.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

####Ploting
Cumulative_Variance <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/12_Cumulative_Variance.txt")
p <- ggplot(Cumulative_Variance, aes(x = PCs, y = Cumulative_Variance)) + geom_line() + geom_point()
p <- p + theme_bw()
p <- p + geom_hline (yintercept = 0.95, linetype="dashed", color = "red")
p <- p + scale_x_continuous(breaks = round(seq(1, 369, by = 20),1))
p <- p + scale_y_continuous(breaks = round (seq (0, 1, by = 0.1), 1))
p <- p + xlab ("PCs") + ylab ("Percentage of Cumulative Variance")
p <- p +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=18, face = "bold"),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs(title="Cumulative Variance Plot")
        
######Loadings Plot
df_x <- as.data.frame(df_PCA$x)
write.table(df_x, file = "/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/11_Data_Analysis_1_PCA_X.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

metabolite_selected <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/14_Data_Analysis_1_PCA_Loading_PC11_PC2.txt")

p <- ggplot(metabolite_selected, aes(x = PC11, y = PC2, label = Metabolites), size=3, vjust="outward")+
  geom_point() + geom_text_repel(max.overlaps = Inf, show.legend  = F)
p <- p + geom_segment(aes(xend=PC11, yend=PC2), x=0, y=0, color="red")
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
