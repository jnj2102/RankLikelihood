%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
%Cholesky Decomposition method
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n50, p25, sparsity = AR4

%set the dimension
p=25;
%set the sample size
n=50;

%these are all true omegas (omega_true) - I will need to loop through the models
[star_model, AR2_model, AR4_model, circle_model, full_model, AR1_model, band_model,...
    fivepercent_model,tenpercent_model, twentypercent_model, twopercent_model,...
    fifteenpercent_model] = model_simulation(p);



%set the true omega
omega_true = AR4_model;
sigma_true = inv(omega_true);

%create the true means. 
mean_true = transpose( linspace(0,5,p)); 

mu = 1; %I'm just giving it a mean of 1 for now but it can be any constant
tau = 1; %I'm giving it a sd of 1 for now but it can be any constant

sigma2 = 1; %I'm making the variance 1 for now but it can be any constant

%set the number of points for the Gauss-Hermite quadrature to 20

N_points = 20;

c = [0;1]; %vector of linear constraints

%save the sigma matrix for R
%save('RankLikelihood_Sigma_true_p25_n50_AR4.mat','sigma_true');

%save the omega matrix for R
%save('RankLikelihood_Omega_true_p25_n50_AR4.mat','omega_true');

%set the number of replications for the simulation
reps = 100;

    y_true_n50_p25 = cell([reps,1]);

    optimalJ_n50_p25_AR4 = cell([reps,1]);
    minK_n50_p25_AR4 = cell([reps,1]);
    Final_AIC_n50_p25_AR4 = cell([reps,1]);
    

 x_matrix_n50_p25 =  cell([reps,1]);
       
F_mat_cell_iters = cell([reps,1]);
g_vec_cell_iters = cell([reps,1]);
A_mat_cell_iters = cell([reps,1]);
LHS_gq_cell_iters = cell([reps,1]);
W_mat_cell_iters = cell([reps,1]);
q_vec_cell_iters = cell([reps,1]);
knot_vector_cell_iters = cell([reps,1]);
index_1_cell_iters = cell([reps,1]);
index_2_cell_iters = cell([reps,1]);
inverse_variance_prior_reduced_cell_iters = cell([reps,1]);
mean_prior_reduced_cell_iters = cell([reps,1]);
Fconstraint_cell_iters = cell([reps,1]);
RHS_gq_cell_iters = cell([reps,1]);
Z_red_cell_iters = cell([reps,1]);
Z_two_cell_iters = cell([reps,1]);

MLE_ald_n50_p25= cell([reps,1]);
    MLE_ev_n50_p25 = cell([reps,1]);
    MLE_stable_n50_p25 = cell([reps,1]);
    MLE_log_n50_p25 = cell([reps,1]);
    initial_value_cell_iters = cell([reps,1]);



%Run the iterations

for iters = 1:reps
   
    fprintf('iters');
    
    rng(100 + iters,'twister'); %set the seed for each replication for reproducibility
      
    % generate the true y using the AR4 model
    y_true = mvnrnd(mean_true,sigma_true,n); 
    
    
    %Save the true y
    y_true_n50_p25{iters} = y_true;
     
   
    [x_matrix, MLE_ald,MLE_ev, MLE_log, MLE_stable] = xfunctions_ranklikelihood(p, y_true);
    
    x_matrix_n50_p25{iters} = x_matrix;
    MLE_ev_n50_p25{iters} = MLE_ev;
    MLE_stable_n50_p25{iters} = MLE_stable;
    MLE_ald_n50_p25{iters} = MLE_ald;
    MLE_log_n50_p25{iters} =  MLE_log;


    
    %save the x matrix as a mat file for R
    
   % save(sprintf('RankLikelihood_Iter_%d_X_matrix_p25_n50_AR4.mat', iters), 'x_matrix');
       
    %First find the optimal J for the replication, then the prior mean and
    %variance and finally the initial values
    
[optimalJ, minK, Final_AIC,F_mat_cell,g_vec_cell,...
Fconstraint_cell,RHS_gq_cell,A_mat_cell,LHS_gq_cell,W_mat_cell,...
q_vec_cell,knot_vector_cell,index_1_cell,index_2_cell,...
inverse_variance_prior_reduced_cell,mean_prior_reduced_cell,...
Z_red_cell,Z_two_cell,initial_value] = PriorBsplines_Initialvalues_corrected(n,p, x_matrix, mu, tau, sigma2,c,...
N_points);

  
    optimalJ_n50_p25_AR4{iters} = optimalJ;
    minK_n50_p25_AR4{iters} = minK;
    Final_AIC_n50_p25_AR4{iters} = Final_AIC;
    
    %Initialize the cells again for the dth predictor

    
    
initial_value_cell_iters{iters} = initial_value;

%I have these iters to run in my matlab function for each iters
F_mat_cell_iters{iters} = F_mat_cell;
g_vec_cell_iters{iters} = g_vec_cell;
Fconstraint_cell_iters{iters} = Fconstraint_cell;
RHS_gq_cell_iters{iters} = RHS_gq_cell;
A_mat_cell_iters{iters} = A_mat_cell;
LHS_gq_cell_iters{iters} = LHS_gq_cell;
W_mat_cell_iters{iters} = W_mat_cell;
q_vec_cell_iters{iters} = q_vec_cell;
knot_vector_cell_iters{iters} = knot_vector_cell;
index_1_cell_iters{iters} = index_1_cell;
index_2_cell_iters{iters} = index_2_cell;
inverse_variance_prior_reduced_cell_iters{iters} = inverse_variance_prior_reduced_cell;
mean_prior_reduced_cell_iters{iters} = mean_prior_reduced_cell;
Z_red_cell_iters{iters} = Z_red_cell;
Z_two_cell_iters{iters}  = Z_two_cell;



end %end of iters loop


save('RankLikelihood_p25_n50_AR4_prior.mat');

    