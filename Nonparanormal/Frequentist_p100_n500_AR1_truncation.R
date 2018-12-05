########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Rank Likelihood
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(huge)
library(orca)
library(pulsar)

source("Frequentist_truncation.R") #call the function



reps <- 100

###########Simulation combination: n=500, p=100, sparsity = AR1#######
result_list <- list()


data <- readMat('RankLikelihood_p100_n500_AR1_prior.mat', package = "R.matlab") #Read in all the data


for (iters in 1:reps) {
  cat('Iteration = ', iters)

set.seed(100 + iters)
  
#Run this in a loop for each iteration 

  
  #read in the Sigma_true 
  
  sigma_true <- data$sigma.true
  

  #read in the omega_true
  
  omega_true <- data$omega.true
  
  #read in the x matrix (observed data)
  
xmat <- data$x.matrix.n500.p100[[iters]][[1]]

p = ncol(xmat)
n = nrow(xmat)



#True edge matrix for the model 
edge_matrix_true <- matrix(0, nrow = p, ncol = p)

for (j in 1:p) {
  for (k in 1:p) {
    
    if ( omega_true[j,k] == 0) {
      edge_matrix_true[j,k] = 0
    } else {
      edge_matrix_true[j,k] = 1
    }
    
  }
}

#find the upper indices
indmx = matrix(1:p^2, nrow = p, ncol = p)

upperind_diag = which(triu(indmx)>0)  #include the diagonal

upperind = which(triu(indmx,1)>0)  #do not include the diagonal



#nonparanormal truncation with graphical lasso
result_list[[iters]] <- Frequentist_truncation(omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind)


}

###################################################################

#Save the data for latex later

save.image(file = "Frequentist_p100_n500_AR1_truncation.rdata")
     
