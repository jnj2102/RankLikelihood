####################################################
# Boxplots for the RankLikelihood paper
# Author: Jami Jackson Mulgrave
#
#
######################################################


#clear the workspace
rm(list = ls())


fullTable <- read.csv(file = 'RankLikelihood_Boxplot_Loss_p25_n50.csv', header = TRUE, sep = ',')

fullTable <- as.data.frame(fullTable)

fullTable$Dimension <- as.factor(fullTable$Dimension)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_tenpercent.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
BLoss_value[i] = result_list[[i]]$bounded_loss
FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision

}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_AR4.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n50_AR1.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR4_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR1_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_tenpercent_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_tenpercent_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR1_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p25_n50_AR4_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the p=50 data


fullTable_p50 <- read.csv(file = 'RankLikelihood_Boxplot_Loss_p50_n100.csv', header = TRUE, sep = ',')

fullTable_p50 <- as.data.frame(fullTable_p50)

fullTable_p50$Dimension <- as.factor(fullTable_p50$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p50)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_fivepercent.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR4.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR1.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR4_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR1_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR1_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR4_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)




#Read in the p=100 data


fullTable_p100 <- read.csv(file = 'RankLikelihood_Boxplot_Loss_p100_n500.csv', header = TRUE, sep = ',')

fullTable_p100 <- as.data.frame(fullTable_p100)

fullTable_p100$Dimension <- as.factor(fullTable_p100$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p100)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_twopercent.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_AR4.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p100_n500_AR1.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, BDGraph_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR4_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR1_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_twopercent_skeptic.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('SKEPTIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_twopercent_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR1_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(1)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(ELoss_value, BLoss_value, FrobLoss_value, Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p100_n500_AR4_truncation.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(4)', number_elements)
Method = rep('Truncation', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

save(fullTable, file = 'RankLikelihood_Boxplot_Loss_table.rdata')
