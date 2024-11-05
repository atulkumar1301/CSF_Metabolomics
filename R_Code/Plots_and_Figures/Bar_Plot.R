#! /Library/Frameworks/R.framework/Versions/4.2/Resources/bin/Rscript
library(data.table)
library (ggplot2)

df <- fread ("Result_Mean_Metabolite_Level_Pathology.txt")

p <- ggplot(data=df, aes(x=Model, y=`t Value`)) +
  geom_bar(stat="identity")

p + scale_x_discrete(limits=c(expression("A"*beta), expression("A"*beta*" adjusted for other co-pathology"), "TauPET", "TauPET adjusted for other co-pathology",
                              expression (alpha*"Syn SAA"), expression (alpha*"Syn SAA adjusted for other co-pathology"), "WML", "TWML adjusted for other co-pathology"))
p
