########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Cholesky Decomposition
#Author: Jami Jackson Mulgrave
#
########################

Frequentist_truncation <- function(omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind) {
  
ptm <- proc.time()

Ymat_npn = huge.npn(xmat, npn.func = "truncation") # Nonparanormal truncation

algorithm_time <- proc.time() - ptm  #I would use the user value for the algorithm time.  


# ### Use G-STARS for Model Selection #####

lambda.min.ratio = 0.01
nlambda = 100

Ymat_npn_scale = scale(Ymat_npn) #scaling is not necessary to find the correlation matrix but it's what is in 
#the huge package function huge.
Ymat_npn_scale_cor = cor(Ymat_npn_scale)

S = Ymat_npn_scale_cor #correlation matrix


lambda.max = max(max(S - diag(p)), -min(S - diag(p)))  #p is the dimension
#lambda.min.ratio:  it is the smallest value for lambda, as a fraction
#of the uppperbound (MAX) of the regularization/thresholding parameter which
#makes all estimates equal to 0.  lambda.max makes all estimates equal to zero
# The program can automatically generate lambda
# as a sequence of length = nlambda starting from MAX to lambda.min.ratio*MAX
# in log scale. The program can automatically generate lambda as a sequence of length
# = nlambda, which makes the sparsity level of the graph path increases from 0
# to lambda.min.ratio evenly.The default value is 0.1 when method = "mb" or
# "glasso"
lambda.min = lambda.min.ratio * lambda.max
lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))



hugeargs <- list(lambda = lambda,
                 method = "glasso", cov.output = TRUE)

ans_stars_gcd_glasso = pulsar(Ymat_npn, fun = huge::huge, fargs = hugeargs, c('stars', 'gcd'),
                              lb.stars=TRUE, ub.stars=TRUE
                   ,rep.num=100) #glasso

#G-StARS Criterion:
# Use the  default approach for selecting the optimal index, based on the gcd+StARS criterion: 
#choose the minimum gcd summary statistic between the upper and lower StARS bounds.
lambda_selected_GStars_index <- get.opt.index(ans_stars_gcd_glasso, 'gcd')

#Refit function has the inverse covariance matrix.  Use the selected lambda to find the optimal inverse covariance
#matrix


fit.ans_stars_gcd_glasso = refit(ans_stars_gcd_glasso) #refits using the glasso function in huge package with the previous criteria.
#using the lambda selected from stars and gcd criterion, called G-StARS.


# #optimal precision and covariance matrix

opt_precisionMat_stars_gcd_glasso= fit.ans_stars_gcd_glasso$est$icov[[lambda_selected_GStars_index]]
opt_covarianceMat_stars_gcd_glasso= fit.ans_stars_gcd_glasso$est$cov[[lambda_selected_GStars_index]]

edge_matrix_est_stars_gcd_glasso =  fit.ans_stars_gcd_glasso$est$path[[lambda_selected_GStars_index]]

####Original StARS criterion
# 
# lambda_selected_stars_gcd_glasso = lambda[ans_stars_gcd_glasso$stars$opt.index]
# 
# lambda_selected_stars_gcd_glasso_index = ans_stars_gcd_glasso$stars$opt.index
# 
# 
# # #optimal precision and covariance matrix
# 
# opt_precisionMat_stars_gcd_glasso= fit.ans_stars_gcd_glasso$est$icov[[lambda_selected_stars_gcd_glasso_index]]
# opt_covarianceMat_stars_gcd_glasso= fit.ans_stars_gcd_glasso$est$cov[[lambda_selected_stars_gcd_glasso_index]]
# 
# edge_matrix_est_stars_gcd_glasso =  fit.ans_stars_gcd_glasso$est$path[[lambda_selected_stars_gcd_glasso_index]]
# 

# 
# #trace is the sum of the diagonal elements. 
# 
# 
# #Use the optimal precision matrix for the entropy loss

entropy_loss  = sum(diag(opt_precisionMat_stars_gcd_glasso%*% sigma_true)) -
  log(det(opt_precisionMat_stars_gcd_glasso%*%sigma_true)) - p;


#Frobenius loss
FrobeniusLoss_precision = (sum(diag(crossprod(opt_precisionMat_stars_gcd_glasso - omega_true)))) 


FrobeniusLoss_covariance = (sum(diag(crossprod(opt_covarianceMat_stars_gcd_glasso - sigma_true)))) 

#Bounded loss
bounded_loss = 1/(p^2) * sum(sum(abs(opt_precisionMat_stars_gcd_glasso - omega_true)))


#Find the TP, TN, FP, and FN
TP_matrix_stars_gcd_glasso= edge_matrix_est_stars_gcd_glasso[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_stars_gcd_glasso= sum(TP_matrix_stars_gcd_glasso) #the colon sums all elements in the matrix


TN_matrix_stars_gcd_glasso= edge_matrix_est_stars_gcd_glasso[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_stars_gcd_glasso= sum(TN_matrix_stars_gcd_glasso) #the colon sums all elements in the matrix


FP_matrix_stars_gcd_glasso= edge_matrix_est_stars_gcd_glasso[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_stars_gcd_glasso= sum(FP_matrix_stars_gcd_glasso) #the colon sums all elements in the matrix


FN_matrix_stars_gcd_glasso= edge_matrix_est_stars_gcd_glasso[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_stars_gcd_glasso = sum(FN_matrix_stars_gcd_glasso) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_stars_gcd_glasso= TN_stars_gcd_glasso/(TN_stars_gcd_glasso + FP_stars_gcd_glasso)

SE_freq_stars_gcd_glasso= TP_stars_gcd_glasso/(TP_stars_gcd_glasso + FN_stars_gcd_glasso)

temp1_stars_gcd_glasso= as.numeric((TP_stars_gcd_glasso + FP_stars_gcd_glasso)*
                                     (TP_stars_gcd_glasso + FN_stars_gcd_glasso))

temp2_stars_gcd_glasso= as.numeric((TN_stars_gcd_glasso + FP_stars_gcd_glasso)*
                                     (TN_stars_gcd_glasso + FN_stars_gcd_glasso))

MCC_freq_stars_gcd_glasso= as.numeric(((TP_stars_gcd_glasso * TN_stars_gcd_glasso) -
                                         (FP_stars_gcd_glasso * FN_stars_gcd_glasso))/
                                        sqrt(temp1_stars_gcd_glasso*temp2_stars_gcd_glasso))



return(list(algorithm_time = algorithm_time,
            TP_stars_gcd_glasso=TP_stars_gcd_glasso,
            TN_stars_gcd_glasso = TN_stars_gcd_glasso, FP_stars_gcd_glasso = FP_stars_gcd_glasso, 
            FN_stars_gcd_glasso = FN_stars_gcd_glasso,
            MCC_freq_stars_gcd_glasso = MCC_freq_stars_gcd_glasso,
            SP_freq_stars_gcd_glasso = SP_freq_stars_gcd_glasso,
            SE_freq_stars_gcd_glasso = SE_freq_stars_gcd_glasso,
            entropy_loss = entropy_loss,
            FrobeniusLoss_precision = FrobeniusLoss_precision,
            FrobeniusLoss_covariance = FrobeniusLoss_covariance,

opt_precisionMat_stars_gcd_glasso= opt_precisionMat_stars_gcd_glasso,
opt_covarianceMat_stars_gcd_glasso= opt_covarianceMat_stars_gcd_glasso,

edge_matrix_est_stars_gcd_glasso =  edge_matrix_est_stars_gcd_glasso,

            ans_stars_gcd_glasso = ans_stars_gcd_glasso,
           fit_ans_stars_gcd_glasso = fit.ans_stars_gcd_glasso,
						lambda_selected_GStars_index = lambda_selected_GStars_index,
						bounded_loss = bounded_loss))

}
