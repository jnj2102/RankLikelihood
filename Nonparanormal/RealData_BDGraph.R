########################
#Code to run the BDMCMC method using the BDGraph package
#For the Simulation section of my paper
#RankLikelihood
#REal Data Example
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(BDgraph)
library(Matrix)


source("BDGraph_copula.R") #call the function




#Data from BDGraph package

data <- readMat('BDGraph_dataset.mat', package = "R.matlab") #Read in all the data




set.seed(5001) #different seed from regression approach.
  

xmat <- data$bdgraphdata


p = dim(xmat)[2]

#Just put identity for the truth - not using those results.
  
  sigma_true <- diag(p)
  

  #read in the omega_true
  
  omega_true <- diag(p)
  
  


#nonparanormal truncation with graphical lasso
result_list <- BDGraph_copula(omega_true, sigma_true, xmat,p)




###################################################################

#Save the data for latex later

save.image(file = "RealData_BDGraph.rdata")
     
writeMat('RealData_BDGraph.mat', edgeMat_BDGraph = result_list$edge_matrix_selected)
