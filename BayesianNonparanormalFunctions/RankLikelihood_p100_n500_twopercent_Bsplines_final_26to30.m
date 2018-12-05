%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=500, p=100, sparsity = twopercent



 c_tune_list = [.1;  1;  10; 100];
 
 ssIters = 26;
 ssIters2 = ssIters;
 
 while ssIters <= 30

      
 for num_iters = ssIters:ssIters2
     
    rng(200 + num_iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', num_iters);
       
for Number_index = 1:4      
%Horseshoe
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_%dto%d_%d.mat', ssIters,ssIters2,Number_index) );

   if c_tune ~= c_tune_list(Number_index)
       fprintf('For p100_n500_twopercent_Bsplines_%dto%d_%d, c_tune is not %d', ssIters,ssIters2,Number_index,c_tune_list(Number_index));
    break;
   else       

    %Now do the B-splines estimation method    
    
    % A_index =1 is already saved.
               
        Omega_Bayes_est =  Omega_Bayes_est_n500_p100_twopercent{A_index, num_iters};
        final_edge_matrix = edge_matrix_n500_p100_twopercent{A_index, num_iters};
        
       mean_Z_Bayes_est = mean_Z_Bayes_est_n500_p100_twopercent{A_index, num_iters};

       
       %%%Find the MLE using a convex optimization problem for eBIC

     %%% MLE with Bglasso
     
[BIC]= FindBIC(Omega_Bayes_est,final_edge_matrix,mean_Z_Bayes_est,n,p);

     BIC_matrix(Number_index, num_iters) = BIC;
     
         %Reconstruct the matrices and such that I plan to search through later
     
     SP_matrix_bglasso_tmp(Number_index,num_iters)  = SP_matrix_bglasso(A_index, num_iters);
   SE_matrix_bglasso_tmp(Number_index,num_iters)  = SE_matrix_bglasso(A_index, num_iters);
   MCC_matrix_bglasso_tmp(Number_index,num_iters)  = MCC_matrix_bglasso(A_index, num_iters);
   total_time_n500_p100_twopercent_tmp(Number_index,num_iters) = total_time_n500_p100_twopercent(A_index, num_iters);
   entropy_loss_n500_p100_twopercent_tmp(Number_index,num_iters) = entropy_loss_n500_p100_twopercent(A_index, num_iters);
  bounded_loss_n500_p100_twopercent_tmp(Number_index,num_iters) = bounded_loss_n500_p100_twopercent(A_index, num_iters);
   Frobenius_norm_precision_n500_p100_twopercent_tmp(Number_index,num_iters) = Frobenius_norm_precision_n500_p100_twopercent(A_index, num_iters);
  Frobenius_norm_covariance_n500_p100_twopercent_tmp(Number_index,num_iters) = Frobenius_norm_covariance_n500_p100_twopercent(A_index, num_iters);
Omega_Bayes_est_n500_p100_twopercent_tmp{Number_index,num_iters} = Omega_Bayes_est_n500_p100_twopercent{A_index, num_iters};
Sigma_Bayes_est_n500_p100_twopercent_tmp{Number_index,num_iters} = Sigma_Bayes_est_n500_p100_twopercent{A_index, num_iters};
edge_matrix_n500_p100_twopercent_tmp{Number_index,num_iters} = edge_matrix_n500_p100_twopercent{A_index, num_iters};
TP_bglasso_n500_p100_twopercent_tmp(Number_index,num_iters) = TP_bglasso_n500_p100_twopercent(A_index, num_iters);
TN_bglasso_n500_p100_twopercent_tmp(Number_index,num_iters) = TN_bglasso_n500_p100_twopercent(A_index, num_iters);
 FP_bglasso_n500_p100_twopercent_tmp(Number_index,num_iters) = FP_bglasso_n500_p100_twopercent(A_index, num_iters);
  FN_bglasso_n500_p100_twopercent_tmp(Number_index,num_iters) =  FN_bglasso_n500_p100_twopercent(A_index, num_iters);
mean_Z_Bayes_est_n500_p100_twopercent_tmp{Number_index,num_iters} = mean_Z_Bayes_est_n500_p100_twopercent{A_index, num_iters};

   
    end
end
end


    
    %Find the minimum BIC and keep the data that pertains to it for further
    %analysis 
for num_iters = ssIters:ssIters2
[minBIC, ~ ] = min(min(BIC_matrix(:,num_iters)));

[rowMinBIC, ~] = find(BIC_matrix(:,num_iters) == minBIC);

%save the data with the minimum BIC for further analysis
 
 SP_matrix_finalanalysis(num_iters)  =  SP_matrix_bglasso_tmp(rowMinBIC,  num_iters);
 SE_matrix_finalanalysis(num_iters) =  SE_matrix_bglasso_tmp(rowMinBIC,  num_iters);
 MCC_matrix_finalanalysis(num_iters) =  MCC_matrix_bglasso_tmp(rowMinBIC,  num_iters);
 
BIC_matrix_finalanalysis(num_iters) =   BIC_matrix(rowMinBIC,  num_iters);

 total_time_finalanalysis(num_iters) =   total_time_n500_p100_twopercent_tmp(rowMinBIC,  num_iters);
entropy_loss_finalanalysis(num_iters)  = entropy_loss_n500_p100_twopercent_tmp(rowMinBIC, num_iters) ;
  bounded_loss_finalanalysis(num_iters)  = bounded_loss_n500_p100_twopercent_tmp(rowMinBIC, num_iters) ;
   Frobenius_norm_precision_finalanalysis(num_iters)  = Frobenius_norm_precision_n500_p100_twopercent_tmp(rowMinBIC, num_iters) ;
  Frobenius_norm_covariance_finalanalysis(num_iters)  = Frobenius_norm_covariance_n500_p100_twopercent_tmp(rowMinBIC, num_iters) ;

 Omega_Bayes_est_finalanalysis{num_iters} =   Omega_Bayes_est_n500_p100_twopercent_tmp{rowMinBIC,  num_iters};
Sigma_Bayes_est_finalanalysis{num_iters} = Sigma_Bayes_est_n500_p100_twopercent_tmp{rowMinBIC,  num_iters};
edge_matrix_finalanalysis{num_iters} =  edge_matrix_n500_p100_twopercent_tmp{rowMinBIC,  num_iters};
   TP_finalanalysis(num_iters) =  TP_bglasso_n500_p100_twopercent_tmp(rowMinBIC,  num_iters);
  TN_finalanalysis(num_iters) =  TN_bglasso_n500_p100_twopercent_tmp(rowMinBIC,  num_iters);
  FP_finalanalysis(num_iters) =  FP_bglasso_n500_p100_twopercent_tmp(rowMinBIC,  num_iters);
   FN_finalanalysis(num_iters) = FN_bglasso_n500_p100_twopercent_tmp(rowMinBIC,  num_iters);
   
  mean_Z_Bayes_est_finalanalysis{num_iters} = mean_Z_Bayes_est_n500_p100_twopercent_tmp{rowMinBIC,  num_iters};
  
end



ssIters = ssIters2 +1;
ssIters2 = ssIters;

 end
 
 
 
%Save one final file with all of the data.
save('RankLikelihood_p100_n500_twopercent_Bsplines_final_26to30.mat', '-v7.3');
