function [SP_bglasso, SE_bglasso, MCC_bglasso,...
   edge_matrix_bglasso, inverse_correlation_Bayes_est,...
    total_time,mean_Z_Bayes_est,Frobenius_norm_precision,...
    Frobenius_norm_covariance,Sigma_Bayes_est, ...
    entropy_loss, TP_bglasso, TN_bglasso, FP_bglasso,FN_bglasso,...
    bounded_loss] = RankLikelihood_CholeskyDecomp_MCMC_horseshoe_mem(n,p,...
     x_matrix, omega_true, sigma_true, capital_A, a0, b0, c_tune)

%Function to run the B-splines approach with the Cholesky decomposition for
%the Bayesian nonparanormal.

%Author: Jami Jackson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indmx = reshape(1:p^2,p,p); 
upperind = indmx(triu(indmx,1)>0);  %does not include the diagonal
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ind_noi_all = zeros(p-1,p);

for i = 1:p
       if i==1  
       ind_noi = (2:p)'; 
      elseif i==p
       ind_noi = (1:p-1)'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
       
end

   
for d = 1:p

    [~,~,ranks]=unique(x_matrix(:,d)); %ascending order
    ranks_matrix(:,d) = ranks;  %this matches the R

end

ranks_matrix_transpose = ranks_matrix';  %this is p x n
   

Z = norminv((ranks_matrix_transpose/(n+1))',0,1);  %transpose again for n x p

%if all eigenvalues of the correlation matrix 
    %isn't positive, use nearestSPD  


if all(eig(corrcoef(Z)) > 0) ~= 1 
    
inverse_correlation_regularized = invChol_mex(nearestSPD(corrcoef(Z)));

else 
  inverse_correlation_regularized = invChol_mex(corrcoef(Z));

end                                 
 
 
%initial values for the Cholesky algorithm

b_offdiag_initial = cell([p-1,1]);


for m = 1:(p-1) 
    
    b_temp = ones([p-m,1]);
    b_offdiag_initial{m} = b_temp;
end
b_offdiag_cell = b_offdiag_initial;
inverse_sigma_squared_vector = ones([p,1]);
inverse_lambda_squared_vector = ones([p-1,1]);


inverse_a_offdiag_vector = inverse_lambda_squared_vector;
inverse_c_offdiag_cell = b_offdiag_cell;


%Start the MCMC

 
   
burnin = 5000;
        

num_samples = 15000; 

      
    complete_Z_Bayes_est_sims = zeros([n,p,num_samples - burnin]);
    complete_inverse_correlation_sims = zeros([p,p,num_samples - burnin]);
     complete_Sig_total_sims = zeros([p,p,num_samples - burnin]);
    
%Step 5 Truncated Multivariate Normal sampling using HMC
tic
for iters_num = 1:num_samples 
    rng(2000 + iters_num,'twister'); %set the seed for each replication for reproducibility

    
    %updating Z directly in the function saves time.
[Z] =  BayesianRankLikelihood(p,ind_noi_all,Z, inverse_correlation_regularized,n,...
    x_matrix);





[inverse_lambda_squared_vector,inverse_a_offdiag_vector,b_offdiag_cell,inverse_c_offdiag_cell,...
    inverse_sigma_squared_vector,Omega_matrix] = CholeskyDecomp_MCMC_horseshoe(Z,p,...
    inverse_lambda_squared_vector, inverse_sigma_squared_vector,...
    b_offdiag_cell,inverse_a_offdiag_vector, capital_A,inverse_c_offdiag_cell, a0, b0,n,...
    c_tune);

Omega_matrix = 1/2*(Omega_matrix + Omega_matrix'); %do this for numerical stability.

%take the inverse to get sigma
Sig = invChol_mex(Omega_matrix);  %this is the easiest way but perhaps I can avoid this inverse.

Diag_sqrt_covariance = diag(sqrt(diag(Sig)));

%now find the inverse correlation matrix

inverse_correlation_regularized = Diag_sqrt_covariance*Omega_matrix*Diag_sqrt_covariance;

%Only save the samples after burnin
if iters_num > burnin
    %create the corresponding index
    
   save_index = iters_num - burnin;

complete_Sig_total_sims(:,:,save_index) = Sig;


complete_inverse_correlation_sims(:,:,save_index) = inverse_correlation_regularized;

complete_Z_Bayes_est_sims(:,:,save_index) = Z; %these cost some time

end



end %end of number of iterations
  total_time=  toc;  %%END OF MCMC %%

   
%Do the burnin for the Z's.  These Z's are Z ~ N(0, Sigma). Not the edge
%matrix.

%complete_Z_Bayes_est_sims = Bayes_est_z_iters(:,:,burnin+1:end);
%Do burnin for omega

%complete_inverse_correlation_sims = inverse_correlation_sims_chain_tmp(:,:, burnin+1:end); %these are the omegas I'm testing 
                                    %for Rhat

% complete_Sig_total_sims = Sig_sims_chain_tmp(:,:,burnin+1:end);

%Summaries
%Bayes estimate of inverse_correlation
inverse_correlation_Bayes_est = mean(complete_inverse_correlation_sims,3);

Sigma_Bayes_est = mean(complete_Sig_total_sims,3);  


%Bayes estimate of Z transformed variables
mean_Z_Bayes_est = mean(complete_Z_Bayes_est_sims,3);

%turn the covariance matrix into a correlation matrix.
    sigma_true_correlation = corrcov(sigma_true); 
    

%Find entropy loss using the inverse correlation matrix since we are
%working with ranks
entropy_loss =  trace(inverse_correlation_Bayes_est*...
            sigma_true_correlation) -...
        log(det(inverse_correlation_Bayes_est*sigma_true_correlation)) - p;

    %Find the inverse of the true correlation matrix
    inverse_sigma_true_correlation = invChol_mex(sigma_true_correlation);
    

    %Find the bounded loss
    
    bounded_loss = 1/(p^2) * sum(sum(abs(inverse_correlation_Bayes_est -...
        inverse_sigma_true_correlation)));
    
    
    %Find the Frobenius Norm 
    Frobenius_norm_precision = trace((inverse_correlation_Bayes_est -...
        inverse_sigma_true_correlation)'*(inverse_correlation_Bayes_est -...
        inverse_sigma_true_correlation));
    
    %make my Sigma_Bayes_est covariance matrix into a correlation matrix
    Sigma_Bayes_est_correlation = corrcov(Sigma_Bayes_est);
    
        Frobenius_norm_covariance = trace((Sigma_Bayes_est_correlation -...
        sigma_true_correlation)'*(Sigma_Bayes_est_correlation -...
        sigma_true_correlation));
    
    
%Now calculate SE, SP, MCC, L1 and give the final estimated edge matrix
                 
%True edge matrix for the model 
edge_matrix_true = zeros(p,p); %1 is edge is present 0 is edge is absent.

for j = 1:p
    for k = 1:p
        if omega_true(j,k) == 0
            edge_matrix_true(j,k) = 0;
        else
            edge_matrix_true(j,k) = 1;
        end
    end
end


%0-1 Loss procedure
partial_corr = zeros(p,p);
 
partial_corr_Cholesky = zeros(p,p);
edge_matrix = zeros(p,p);
 
 
[~,~,num_replicates] = size(complete_inverse_correlation_sims);
 
complete_edge_matrix_total_sims = zeros([p,p, num_replicates]);
 

  %standard conjugate Wishart distribution(3, I_p)
nu = 3;
diag_wishart = eye(p);

nu_post = n + nu;  %degrees of freedom of the posterior


for replicates = 1:num_replicates
   
    Z_temp = complete_Z_Bayes_est_sims(:,:,replicates);
   
    S_temp = Z_temp'*Z_temp;
   
scale_post = invChol_mex(diag_wishart + S_temp); %scale matrix of the posterior
 
 
omega_post = scale_post*nu_post; %mean of posterior of Wishart.  I can't
% the backsplace \ because nu_post is a scalar
 
 
for j = 1:p
    for k = 1:p
        partial_corr(j,k) = -omega_post(j,k)./sqrt(omega_post(j,j)*omega_post(k,k));
    end
end
 
 
Omega_matrix_temp = complete_inverse_correlation_sims(:,:,replicates);
   
for j = 1:p
    for k = 1:p
        partial_corr_Cholesky(j,k) = -Omega_matrix_temp(j,k)./sqrt(Omega_matrix_temp(j,j)*Omega_matrix_temp(k,k));
    end
end
 
pi_tilde = partial_corr_Cholesky ./partial_corr;
                                                                   
for j = 1:p
    for k = 1:p
        if pi_tilde(j,k) > 0.5
            edge_matrix (j,k) = 1;
        else
            edge_matrix (j,k) = 0;
        end
    end
end
 
complete_edge_matrix_total_sims(:,:,replicates) = edge_matrix;
 
end
 
 
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication

%Median probability model
 edge_matrix_bglasso = mean(complete_edge_matrix_total_sims,3)>.5;  %put .5 here because half the time we get 1, have the time 0


  
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication

TP_matrix_bglasso = edge_matrix_bglasso(upperind) == 1 & edge_matrix_true(upperind) == 1;
TP_bglasso = sum(TP_matrix_bglasso(:)); %the colon sums all elements in the matrix


TN_matrix_bglasso = edge_matrix_bglasso(upperind) == 0 & edge_matrix_true(upperind) == 0;
TN_bglasso = sum(TN_matrix_bglasso(:)); 

FP_matrix_bglasso = edge_matrix_bglasso(upperind) == 1 & edge_matrix_true(upperind) == 0;
FP_bglasso = sum(FP_matrix_bglasso(:)); 

FN_matrix_bglasso = edge_matrix_bglasso(upperind) == 0 & edge_matrix_true(upperind) == 1;
FN_bglasso = sum(FN_matrix_bglasso(:)); 
%Find Specificity, Sensitivity, Matthews Correlation Coefficient 

SP_bglasso = TN_bglasso/(TN_bglasso + FP_bglasso);

SE_bglasso = TP_bglasso/(TP_bglasso + FN_bglasso);

MCC_bglasso = ((TP_bglasso*TN_bglasso) - (FP_bglasso*FN_bglasso))/sqrt((TP_bglasso+FP_bglasso)*(TP_bglasso+FN_bglasso)*(TN_bglasso+FP_bglasso)*(TN_bglasso+FN_bglasso));
  



end

