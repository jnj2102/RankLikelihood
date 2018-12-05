####################################################
# Boxplots for the RankLikelihood paper
# Author: Jami Jackson Mulgrave
#
#
######################################################


#clear the workspace
rm(list = ls())


fullTable <- read.csv(file = 'RankLikelihood_Boxplot_Edges_p25_n50.csv', header = TRUE, sep = ',')

fullTable <- as.data.frame(fullTable)

fullTable$Dimension <- as.factor(fullTable$Dimension)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_tenpercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
SP_value[i] = result_list[[i]]$SP
SE_value[i] = result_list[[i]]$SE
MCC_value[i] = result_list[[i]]$MCC

}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_AR1.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_AR4.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR4_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR4_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR1_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR1_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_tenpercent_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_tenpercent_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the p=50 data


fullTable_p50 <- read.csv(file = 'RankLikelihood_Boxplot_Edges_p50_n100.csv', header = TRUE, sep = ',')

fullTable_p50 <- as.data.frame(fullTable_p50)

fullTable_p50$Dimension <- as.factor(fullTable_p50$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p50)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_fivepercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR1.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR4.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR4_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR4_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR1_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR1_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the p=100 data


fullTable_p100 <- read.csv(file = 'RankLikelihood_Boxplot_Edges_p100_n500.csv', header = TRUE, sep = ',')

fullTable_p100 <- as.data.frame(fullTable_p100)

fullTable_p100$Dimension <- as.factor(fullTable_p100$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p100)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_twopercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_AR1.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR1', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_AR4.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR4_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR4_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR4', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR1_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR1_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_twopercent_skeptic.rdata')

#skeptic model

SP_value_skeptic <- c()
SE_value_skeptic <- c()
MCC_value_skeptic <- c()

for (i in 1: reps) {
  SP_value_skeptic[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_skeptic[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_skeptic[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_skeptic, SE_value_skeptic, MCC_value_skeptic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

load('Frequentist_p100_n500_twopercent_truncation.rdata')

#truncation model

SP_value_truncation <- c()
SE_value_truncation <- c()
MCC_value_truncation <- c()

for (i in 1: reps) {
  SP_value_truncation[i] = result_list[[i]]$SP_freq_stars_gcd_glasso
  SE_value_truncation[i] = result_list[[i]]$SE_freq_stars_gcd_glasso
  MCC_value_truncation[i] = result_list[[i]]$MCC_freq_stars_gcd_glasso
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_truncation, SE_value_truncation, MCC_value_truncation)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

save(fullTable, file = 'RankLikelihood_Boxplot_Edges_table.rdata')
