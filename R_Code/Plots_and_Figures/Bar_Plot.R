#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library (ggplot2)

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
  theme(axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=14),
        axis.text.y = element_text(family = "serif", size=14)) +
  labs(tag = "*adjusted for other co-pathology") +
  theme(plot.tag.position = c(0.15, 0.02))
p2
