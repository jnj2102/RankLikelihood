function [Y_updated_matrix] =  BayesianRankLikelihood(p,ind_noi_all,Y_initial, omega_initial,n,...
    x_matrix)

%Function to find the Bayesian rank likelihood
%Author: Jami Jackson 
%Input:
%Output:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for d= 1:p
       
   ind_noi_index = ind_noi_all(:,d);
   
  variance_likelihood = 1/omega_initial(d,d);
   
   omega_removed =  omega_initial(d, ind_noi_index);
   
  [~, sort_index] = sort(x_matrix(:,d)); 
  
mean_likelihood = (-variance_likelihood*omega_removed)*...
        (Y_initial(:,ind_noi_index)'); %this doesn't include the dth column
 

    for i = 1:n
         r_index = sort_index(i);  %current order statistic index or position
         
      % Add the highest value of (+infinity)
       if   i == n
            Y_high = Inf;
       else
            Y_high = Y_initial(sort_index(i + 1) , d);  
       end
   
          %%%   Add the lowest value (-infinity)
       
       if  i == 1
           Y_low = -Inf;
       else
           Y_low = Y_initial(sort_index(i - 1),d); 
       end

    
    % If you wish to simulate a random variable
% 'Z' from the non-standard Gaussian N(m,s^2)
% conditional on l<Z<u, then first simulate
% X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;

truncated_temp = trandn((Y_low-mean_likelihood(r_index))/sqrt(variance_likelihood),...
    (Y_high-mean_likelihood(r_index))/sqrt(variance_likelihood));

Y_new_element = mean_likelihood(r_index) + sqrt(variance_likelihood)*truncated_temp;
    
         Y_initial(r_index,d) = Y_new_element;  %now do the next sample with this newly 
                                                %sampled Y element.

    end

    
end

Y_updated_matrix = Y_initial;  %All of the elements are now updated.






















end