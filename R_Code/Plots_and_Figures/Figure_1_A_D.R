library(data.table)
library(ggpubr)
library(ggplot2)
library(magick)
library(pdftools)
library(cowplot)

##############
Cumulative_Variance <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_1/PCA/12_Cumulative_Variance.txt")
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
        axis.text.x = element_text(family = "serif", size=14, angle = 90),
        axis.text.y = element_text(family = "serif", size=14),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank()) + labs (title = "A")


#########
pdf_img <- pdf_render_page ("~/OneDrive - University of Eastern Finland/Work/Projects/CSF_Metabolomics/Manuscript/New_Figures/metabolomics.pdf", dpi = 300)
img <- image_read (pdf_img)
img_plot <- ggdraw () + draw_image(img)

#####

df <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/mQTL_pQTL/Proetin_Vs_Metabolite_Dynamics.txt")
p1 <- ggscatter(df, x = "mean_standardized_protein_level", y = "mean_standardized_metabolomic_level", color = "red",
                size = 2, alpha = 0.6, ggtheme = theme_bw(), add = "reg.line", conf.int = TRUE, 
                cor.method = "spearman", add.params = list(color = "black", fill = "lightgray")) +
  stat_cor(method = "spearman", r.accuracy = 0.1)
p1 <- p1 + labs(title = "C", x=expression("Mean Standardized Proteome Level"), y=expression("Mean Standardized Metabolome Level"))
p1 <- p1 + scale_x_continuous(breaks = round(seq(-1, 3, by = 0.5),1))
p1 <- p1 + scale_y_continuous(breaks = round (seq (-1, 5, by = 0.5), 1))
p1 <- p1 +
  theme(plot.title = element_text(family = "serif", size=18, face = "bold"),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=14),
        axis.text.y = element_text(family = "serif", size=14),
        panel.background = element_blank())

##################

df <- fread ("~/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Work/Projects/CSF_Metabolomics/Analyses_2/Mean_Metabolite_Level_Pathology/Result_Mean_Metabolite_Level_Pathology.txt")

p2 <- ggplot(data=df, aes(x=Model, y=`t Value`)) +
  geom_bar(stat="identity", width = 0.5)

p2 <- p2 + theme_minimal()


p2 <- p2 + scale_x_discrete(breaks=c("Abeta","Abeta Adjusted","Tau PET", "Tau PET Adjusted", "SAA", "SAA Adjusted", "WML", "WML Adjusted",
                                     "Age", "Age Adjusted", "Gender", "Gender Adjusted", "CU_CI", "CU_CI Adjusted"),
                            labels=c(expression ("A"*beta), expression ("A"*beta*"*"), "TauPET", "TauPET*",
                                     expression (alpha*"Syn SAA"), expression(alpha*"Syn SAA*"), "WML", "WML*", "Age", "Age*", "Gender", "Gender*", "CU/CI Status", "CU/CI status*"))
p2 <- p2 + coord_flip()

p2 <- p2 +
  theme(plot.title = element_text(family = "serif", size=18, face = "bold"),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=14),
        axis.text.y = element_text(family = "serif", size=14)) +
  labs(tag = "*adjusted for other co-pathology", title = "D") +
  theme(plot.tag.position = c(0.15, 0.02))

plot <- ggarrange(p, img_plot, p1, p2)
plot
