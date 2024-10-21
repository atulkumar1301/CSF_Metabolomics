# Library
library(fmsb)
library(data.table)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Radar_Plot_Data_T_Value.txt", check.names = FALSE)
df <- data.frame(df, row.names = 1, check.names = FALSE)
op <- par(family = "serif", font=1)
create_beautiful_radarchart <- function(data, color = c("#009E73", "#E7B800", "#FC4E07", "#56B4E9"), 
                                        vlabels = colnames(df), vlcex = 1,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    df, axistype = 1, seg = 8,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.3), plwd = 1, plty = 5,
    # Customize the grid
    cglcol = "black", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
create_beautiful_radarchart(df, caxislabels = c(-6, -4, -2, 0, 2, 4, 6, 8))

legend(
  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#009E73", "#E7B800", "#FC4E07", "#56B4E9"),
  text.col = "black", cex = 1.5, pt.cex = 2
)
