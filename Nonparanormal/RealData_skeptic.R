########################
#Code to run the Frequentist Nonparanormal skeptic
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


source("Frequentist_skeptic.R") #call the function




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

#nonparanormal skeptic with graphical lasso
result_list <- Frequentist_skeptic(omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind)




###################################################################

#Save the data for latex later

save.image(file = "RealData_skeptic.rdata")
#didn't use the function because skeptic estimates all zeros, so there is an error.

#Running Frequentist_skeptic.R manually, I can get the optimal index via;
# lambda_selected_GStars_index <- get.opt.index(result_list$ans_stars_gcd_glasso, 'gcd')
# 
# ans_stars_gcd_glasso$stars$merge[[lambda_selected_GStars_index]]  #Is this the matrix? It might not be the right one
#      
# writeMat('RealData_skeptic.mat', edgeMat_glasso = result_listedge_matrix_est_stars_gcd_glasso)
# 
