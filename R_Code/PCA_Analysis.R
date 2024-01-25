library(data.table)
library(ggfortify)
library(factoextra)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_1/5_Data_Analysis_1_Full_Imputed.txt")
df_PCA <- prcomp (df [,c (30:398)], center = TRUE, scale. = TRUE)
fviz_eig(df_PCA)

fviz_pca_ind(df_PCA,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(df_PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
