%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the RankLikelihood paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%fivepercent p50 n100 ranks 

load('RankLikelihood_p50_n100_fivepercent_ranks_final.mat');

%Method = {'Ranks', 'B-splines', 'BDGraph', 'Truncation', 'SKEPTIC'};
%Loss_Type = {'ELoss', 'BLoss', 'FrobLoss'};
%Sparsity = {'AR4', 'Star', 'AR1', 'Percent'};

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_fivepercent_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [table_p50_n100_fivepercent_ranks];


clearvars -except combine_tables


%AR4 p50 n100 ranks 

load('RankLikelihood_p50_n100_AR4_ranks_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR4'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR4_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR4_ranks];


clearvars -except combine_tables

%AR1 p50 n100 ranks 

load('RankLikelihood_p50_n100_AR1_ranks_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR1'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR1_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR1_ranks];


clearvars -except combine_tables


%fivepercent p50 n100 Bsplines 

load('RankLikelihood_p50_n100_fivepercent_Bsplines_final.mat');


ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_fivepercent_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_fivepercent_Bsplines];


clearvars -except combine_tables


%AR4 p50 n100 Bsplines 

load('RankLikelihood_p50_n100_AR4_Bsplines_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR4'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR4_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR4_Bsplines];


clearvars -except combine_tables

%AR1 p50 n100 Bsplines 

load('RankLikelihood_p50_n100_AR1_Bsplines_final.mat');

ELoss_value = entropy_loss_finalanalysis';
ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_value = bounded_loss_finalanalysis';
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_value = Frobenius_norm_precision_finalanalysis';
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'AR1'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p50_n100_AR1_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p50_n100_AR1_Bsplines];



%Try writing table as a csv file to read into R

writetable(combine_tables,'RankLikelihood_Boxplot_Loss_p50_n100.csv') 

