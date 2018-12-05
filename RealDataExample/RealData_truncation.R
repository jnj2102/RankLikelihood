########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#RankLikelihood
#Real Data Example
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




#Data from BDGraph package

data <- readMat('BDGraph_dataset.mat', package = "R.matlab") #Read in all the data


set.seed(5000)


xmat <- data$bdgraphdata

p = ncol(xmat)
n = nrow(xmat)

#Just put identity for the truth - not using those results.

sigma_true <- diag(p)


#read in the omega_true

omega_true <- diag(p)

#find the upper indices
indmx = matrix(1:p^2, nrow = p, ncol = p)

upperind_diag = which(triu(indmx)>0)  #include the diagonal

upperind = which(triu(indmx,1)>0)  #do not include the diagonal

edge_matrix_true = diag(p)

#nonparanormal truncation with graphical lasso
result_list <- Frequentist_truncation(omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind)




###################################################################

#Save the data for latex later

save.image(file = "RealData_truncation.rdata")
# 
# 
# lambda_selected_stars_gcd_glasso_index = result_list$ans_stars_gcd_glasso$stars$opt.index
# 
# 
# # #optimal precision and covariance matrix
# 
# edge_matrix_est_stars_gcd_glasso =  result_list$fit_ans_stars_gcd_glasso$est$path[[lambda_selected_stars_gcd_glasso_index]]

     
writeMat('RealData_truncation.mat', edgeMat_glasso = result_list$edge_matrix_est_stars_gcd_glasso)

