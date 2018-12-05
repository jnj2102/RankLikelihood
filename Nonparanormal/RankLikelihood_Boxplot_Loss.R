#Script that takes an Rdata file and constructs a facet grid 
##
# Author: Jami Jackson Mulgrave
#################################################



#clear the workspace
rm(list = ls())

#import libraries

library(ggplot2)
# library(gtable)
# library(grid)
# library(magrittr) # for the %>% that I love so well


load('RankLikelihood_Boxplot_Loss_table.rdata')

#Rename the rows from ranks and Bsplines

fullTable$Loss_Type <- gsub(fullTable$Loss_Type, pattern = 'BLoss', replacement = 'Bounded Loss')
fullTable$Loss_Type <- gsub(fullTable$Loss_Type, pattern = 'FrobLoss', replacement = 'Frobenius Loss')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'AR1', replacement = 'AR(1)')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'AR4', replacement = 'AR(4)')
fullTable$Method <- gsub(fullTable$Method, pattern = 'BDGraph', replacement = 'Bayesian Copula')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'tenpercent', replacement = 'Percent')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'fivepercent', replacement = 'Percent')
fullTable$Sparsity <- gsub(fullTable$Sparsity, pattern = 'twopercent', replacement = 'Percent')

fullTable$Method <- gsub(fullTable$Method, pattern = 'truncation', replacement = 'Truncation')
fullTable$Method <- gsub(fullTable$Method, pattern = 'Ranks', replacement = 'Rank Likelihood')
fullTable$Method <- gsub(fullTable$Method, pattern = 'Bsplines', replacement = 'B-splines')

#Remove Frobenius Loss column

#double check there are only the three losses
table(fullTable$Loss_Type)
table(fullTable$Method)


fullTable_subset = fullTable[which(fullTable$Loss_Type == 'Bounded Loss'), ]


#Specify the order of the methods
fullTable_subset$Method <- factor(fullTable_subset$Method , levels = c("Rank Likelihood","B-splines", "Bayesian Copula",
                                                         "SKEPTIC",  "Truncation"))
# 
# Complete_fullTable_subset = subset(Complete_fullTable, Complete_fullTable$Method == 'Bernoulli-Gaussian'|  
#                 Complete_fullTable$Method == 'Horseshoe'| Complete_fullTable$Method ==  'Variational Bayes')
# 
# 
# table(Complete_fullTable_subset$Method)
# table(Complete_fullTable_subset$Loss_Type)
# table(Complete_fullTable_subset$Sparsity)
# 
# Complete_fullTable_subset$Method <- factor(Complete_fullTable_subset$Method , levels = c("Bernoulli-Gaussian","Horseshoe","Variational Bayes","EBIC",
#                                                          "RIC", "StARS","Bayesian Copula"))

bp <- ggplot(fullTable_subset, aes(x=Dimension, y=Loss_value,group = Method)) + 
  geom_boxplot(aes(fill=Method)) #this makes 1 boxplot for 1 method

#I could add an option to make the boxplot narrow: width = 0.2)
#Make the x axis as a factor since it's just 1 dimension (p)

#create a facet grid with 2 variables

#bp + facet_grid(Sparsity ~ Loss_Type) + facet_wrap(~Dimension)

#bp + facet_grid(Sparsity ~ Loss_Type) 

#bp +  facet_wrap(~Dimension,  scale="free")

bp + facet_grid(Sparsity ~  Dimension, scale = "free")  +
  labs(y = expression(paste("Scaled ", L[1], " Loss", sep = " ")))+ theme(strip.text.x = element_blank())

ggsave("RankLikelihood_Boxplot_Loss.jpg")

#bp +  facet_grid(Sparsity ~ Dimension, scale="free")

#+ theme(strip.text.x = element_blank())
#bp +  facet_wrap(Sparsity~Dimension, scale="free")


