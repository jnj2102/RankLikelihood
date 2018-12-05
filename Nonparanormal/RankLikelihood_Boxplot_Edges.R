#Script that takes an Rdata file and constructs a facet grid 
##
# Author: Jami Jackson Mulgrave
#################################################


#clear the workspace
rm(list = ls())

#import libraries

library(ggplot2)


load('RankLikelihood_Boxplot_Edges_table.rdata')
#Rename the rows from ranks and Bsplines


fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'AR1', replacement = 'AR(1)')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'AR4', replacement = 'AR(4)')
fullTable$Method <- gsub(fullTable$Method, pattern = 'BDGraph', replacement = 'Bayesian Copula')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'fivepercent', replacement = 'Percent')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'tenpercent', replacement = 'Percent')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'twopercent', replacement = 'Percent')
fullTable$Method <- gsub(fullTable$Method, pattern = 'truncation', replacement = 'Truncation')
fullTable$Method <- gsub(fullTable$Method, pattern = 'Ranks', replacement = 'Rank Likelihood')
fullTable$Method <- gsub(fullTable$Method, pattern = 'Bsplines', replacement = 'B-splines')




fullTable$Dimension <- as.factor(fullTable$Dimension)

table(fullTable$Sparsity)
table(fullTable$Dimension)

table(fullTable$Method)

#Specify the order of the methods
fullTable$Method <- factor(fullTable$Method , levels = c("Rank Likelihood","B-splines", "Bayesian Copula",
                                                                "SKEPTIC",  "Truncation"))


#Look at MCC

fullTable_MCC <- subset(fullTable, fullTable$Edges_Type == 'MCC')

bp <- ggplot(fullTable_MCC, aes(x=Dimension, y=Edges_value,group = Method)) + 
  geom_boxplot(aes(fill=Method)) #this makes 1 boxplot for 1 method

#I could add an option to make the boxplot narrow: width = 0.2)
#Make the x axis as a factor since it's just 1 dimension (p)

#create a facet grid with 2 variables

#bp + facet_grid(Sparsity ~ Edges_Type) + facet_wrap(~Dimension)

#bp + facet_grid(Sparsity ~ Edges_Type) 

#bp +  facet_wrap(~Dimension,  scale="free")
#y_row <- as_labeller(c('MCC' = "MCC", 'MCC' = "MCC"))
#p + facet_grid(cyl ~ am, labeller = labeller(am = to_string))

bp + facet_grid(Sparsity ~ Dimension, scale = "free") +
  labs(y = 'Matthews Correlation Coefficient') + theme(strip.text.x = element_blank())
#Note that missing is removed. Get an error 
#> sum(is.na(fullTable_MCC$Edges_value))
#[1] 550

ggsave("RankLikelihood_Boxplot_Edges_MCC.jpg")


#Look at SE

fullTable_SE <- subset(fullTable, fullTable$Edges_Type == 'Sensitivity')

bp_SE <- ggplot(fullTable_SE, aes(x=Dimension, y=Edges_value,group = Method)) + 
  geom_boxplot(aes(fill=Method)) #this makes 1 boxplot for 1 method

bp_SE + facet_grid(Sparsity ~ Dimension, scale = "free") +
  labs(y = 'Sensitivity') + theme(strip.text.x = element_blank())
#Note that missing is removed.

ggsave("RankLikelihood_Boxplot_Edges_SE.jpg")


#Look at SP

fullTable_SP <- subset(fullTable, fullTable$Edges_Type == 'Specificity')

bp_SP <- ggplot(fullTable_SP, aes(x=Dimension, y=Edges_value,group = Method)) + 
  geom_boxplot(aes(fill=Method)) #this makes 1 boxplot for 1 method

bp_SP + facet_grid(Sparsity ~ Dimension, scale = "free") +
  labs(y = 'Specificity') + theme(strip.text.x = element_blank())


ggsave("RankLikelihood_Boxplot_Edges_SP.jpg")
#bp +  facet_grid(Sparsity ~ Dimension, scale="free")
#+ theme(strip.text.x = element_blank())
#, labeller = labeller(Sparsity = c("MCC", "MCC", "MCC", "Sensitivity", "Sensitivity", "Sensitivity", "Specificity", "Specificity", "Specificity"))

#bp +  facet_wrap(Sparsity~Dimension, scale="free")
#, labeller = labeller(Edges_Type = y_row)

