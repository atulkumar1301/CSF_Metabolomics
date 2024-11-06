#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library (ggplot2)

df <- fread ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Mean_Metabolite_Level_Pathology/Result_Mean_Metabolite_Level_Pathology.txt")

p <- ggplot(data=df, aes(x=Model, y=`t Value`)) +
  geom_bar(stat="identity", width = 0.5) 
  

p <- p + scale_x_discrete(breaks=c("Abeta","Abeta Adjusted","Tau PET", "Tau PET Adjusted", "SAA", "SAA Adjusted", "WML", "WML Adjusted"),
                     labels=c(expression ("A"*beta), expression ("A"*beta*"*"), "TauPET", "TauPET*",
                              expression (alpha*"Syn SAA"), expression(alpha*"Syn SAA*"), "WML", "WML*"))

p <- p +
  theme(axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=12, angle = 10),
        axis.text.y = element_text(family = "serif", size=12)) +
  labs(tag = "*adjusted for other co-pathology") +
  theme(plot.tag.position = c(0.15, 0.02))
p
