library(data.table)
library(ggplot2)
library(ggpubr)
df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Imputed_Values/Observed_Vs_Missing_Values.txt")
p1 <- ggplot(df, aes(x=A_100001256, fill=Type_A_100001256, color=Type_A_100001256)) +
  geom_histogram(position="dodge")
p1 <- p1 + theme_light()
p1 <- p1 + scale_x_continuous(breaks = round(seq(0, 5, by = 1),1))
p1 <- p1 + scale_y_continuous(breaks = round (seq (0, 180, by = 20), 1))
p1 <- p1 + xlab ("N-acetylphenylalanine") + ylab ("Count")
p1 <- p1 + guides(fill=guide_legend(title="Count"))
p1 <- p1 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="A) N-acetylphenylalanine", face = "bold")
########

p2 <- ggplot(df, aes(x=A_356, fill=Type_A_356, color=Type_A_356)) +
  geom_histogram(position="dodge")
p2 <- p2 + theme_light()
p2 <- p2 + scale_x_continuous(breaks = round(seq(0, 8, by = 1),1))
p2 <- p2 + scale_y_continuous(breaks = round (seq (0, 280, by = 40), 1))
p2 <- p2 + xlab ("Cortisol") + ylab ("Count")
p2 <- p2 + guides(fill=guide_legend(title="Count"))
p2 <- p2 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="B) Cortisol", face = "bold")
#########

p3 <- ggplot(df, aes(x=A_100002500, fill=Type_A_100002500, color=Type_A_100002500)) +
  geom_histogram(position="dodge")
p3 <- p3 + theme_light()
p3 <- p3 + scale_x_continuous(breaks = round(seq(0, 200, by = 20),1))
p3 <- p3 + scale_y_continuous(breaks = round (seq (0, 850, by = 100), 1))
p3 <- p3 + xlab ("Formiminoglutamate") + ylab ("Count")
p3 <- p3 + guides(fill=guide_legend(title="Count"))
p3 <- p3 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="C) Formiminoglutamate", face = "bold")
#########

p4 <- ggplot(df, aes(x=A_825, fill=Type_A_825, color=Type_A_825)) +
  geom_histogram(position="dodge")
p4 <- p4 + theme_light()
p4 <- p4 + scale_x_continuous(breaks = round(seq(0, 450, by = 60),1))
p4 <- p4 + scale_y_continuous(breaks = round (seq (0, 900, by = 100), 1))
p4 <- p4 + xlab ("Uracil") + ylab ("Count")
p4 <- p4 + guides(fill=guide_legend(title="Count"))
p4 <- p4 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="D) Uracil", face = "bold")
#########

p5 <- ggplot(df, aes(x=A_100001175, fill=Type_A_100001175, color=Type_A_100001175)) +
  geom_histogram(position="dodge")
p5 <- p5 + theme_light()
p5 <- p5 + scale_x_continuous(breaks = round(seq(0, 22, by = 3),1))
p5 <- p5 + scale_y_continuous(breaks = round (seq (0, 500, by = 50), 1))
p5 <- p5 + xlab ("1,3-propanediol") + ylab ("Count")
p5 <- p5 + guides(fill=guide_legend(title="Count"))
p5 <- p5 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="E) 1,3-propanediol", face = "bold")
#########

p6 <- ggplot(df, aes(x=A_100020417, fill=Type_A_100020417, color=Type_A_100020417)) +
  geom_histogram(position="dodge")
p6 <- p6 + theme_light()
p6 <- p6 + scale_x_continuous(breaks = round(seq(0, 7, by = 1),1))
p6 <- p6 + scale_y_continuous(breaks = round (seq (0, 150, by = 15), 1))
p6 <- p6 + xlab ("N6,N6-dimethyllysine") + ylab ("Count")
p6 <- p6 + guides(fill=guide_legend(title="Count"))
p6 <- p6 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="F) N6,N6-dimethyllysine", face = "bold")
#########

p7 <- ggplot(df, aes(x=A_339, fill=Type_A_339, color=Type_A_339)) +
  geom_histogram(position="dodge")
p7 <- p7 + theme_light()
p7 <- p7 + scale_x_continuous(breaks = round(seq(0, 16, by = 2),1))
p7 <- p7 + scale_y_continuous(breaks = round (seq (0, 450, by = 50), 1))
p7 <- p7 + xlab ("Glutarate (C5-DC)") + ylab ("Count")
p7 <- p7 + guides(fill=guide_legend(title="Count"))
p7 <- p7 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="G) Glutarate (C5-DC)", face = "bold")
#########

p8 <- ggplot(df, aes(x=A_302, fill=Type_A_302, color=Type_A_302)) +
  geom_histogram(position="dodge")
p8 <- p8 + theme_light()
p8 <- p8 + scale_x_continuous(breaks = round(seq(0, 28, by = 4),1))
p8 <- p8 + scale_y_continuous(breaks = round (seq (0, 450, by = 50), 1))
p8 <- p8 + xlab ("Deoxycholate") + ylab ("Count")
p8 <- p8 + guides(fill=guide_legend(title="Count"))
p8 <- p8 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="H) Deoxycholate", face = "bold")
#########

p9 <- ggplot(df, aes(x=A_999915245, fill=Type_A_999915245, color=Type_A_999915245)) +
  geom_histogram(position="dodge")
p9 <- p9 + theme_light()
p9 <- p9 + scale_x_continuous(breaks = round(seq(0, 5, by = 1),1))
p9 <- p9 + scale_y_continuous(breaks = round (seq (0, 125, by = 25), 1))
p9 <- p9 + xlab ("X-15245") + ylab ("Count")
p9 <- p9 + guides(fill=guide_legend(title="Count"))
p9 <- p9 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="I) X-15245", face = "bold")
#########

p10 <- ggplot(df, aes(x=A_100000096, fill=Type_A_100000096, color=Type_A_100000096)) +
  geom_histogram(position="dodge")
p10 <- p10 + theme_light()
p10 <- p10 + scale_x_continuous(breaks = round(seq(0, 44, by = 5),1))
p10 <- p10 + scale_y_continuous(breaks = round (seq (0, 600, by = 60), 1))
p10 <- p10 + xlab ("4-guanidinobutanoate") + ylab ("Count")
p10 <- p10 + guides(fill=guide_legend(title="Count"))
p10 <- p10 +
  theme(legend.position="left",
        plot.title = element_text(family = "serif", size=12, face = "bold"),
        axis.title.x = element_text(family = "serif", size=12),
        axis.title.y = element_text(family = "serif", size=12),
        axis.text.x = element_text(family = "serif", size=10),
        axis.text.y = element_text(family = "serif", size=10),
        legend.title = element_text(family = "serif", size=14),
        legend.text = element_text(family = "serif", size=14),
        panel.background = element_blank()) + labs(title="J) 4-guanidinobutanoate", face = "bold")


ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, common.legend = TRUE)
