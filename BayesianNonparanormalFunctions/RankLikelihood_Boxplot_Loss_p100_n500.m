%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the RankLikelihood paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%twopercent p100 n500 ranks 

load('RankLikelihood_p100_n500_twopercent_ranks_final.mat');

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

table_p100_n500_twopercent_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [table_p100_n500_twopercent_ranks];


clearvars -except combine_tables


%AR4 p100 n500 ranks 

load('RankLikelihood_p100_n500_AR4_ranks_final.mat');

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

table_p100_n500_AR4_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR4_ranks];


clearvars -except combine_tables

%AR1 p100 n500 ranks 

load('RankLikelihood_p100_n500_AR1_ranks_final.mat');

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

table_p100_n500_AR1_ranks = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR1_ranks];


clearvars -except combine_tables

%twopercent p100 n500 Bsplines 
ssIters = 1;
ssIters2 = 25;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 1;
ssIters2 = 25;

ELoss_value_tmp = entropy_loss_finalanalysis(:);
ELoss_value_1to25 = ELoss_value_tmp(ssIters:ssIters2);
BLoss_value_tmp = bounded_loss_finalanalysis(:);
BLoss_value_1to25 = BLoss_value_tmp(ssIters:ssIters2);
FrobLoss_value_tmp = Frobenius_norm_precision_finalanalysis(:);
FrobLoss_value_1to25 = FrobLoss_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2

ssIters = 26;
ssIters2 = 30;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 26;
ssIters2 = 30;

ELoss_value_tmp = entropy_loss_finalanalysis(:);
ELoss_value_26to30 = ELoss_value_tmp(ssIters:ssIters2);
BLoss_value_tmp = bounded_loss_finalanalysis(:);
BLoss_value_26to30 = BLoss_value_tmp(ssIters:ssIters2);
FrobLoss_value_tmp = Frobenius_norm_precision_finalanalysis(:);
FrobLoss_value_26to30 = FrobLoss_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2

ssIters = 31;
ssIters2 = 75;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 31;
ssIters2 = 75;

ELoss_value_tmp = entropy_loss_finalanalysis(:);
ELoss_value_31to75 = ELoss_value_tmp(ssIters:ssIters2);
BLoss_value_tmp = bounded_loss_finalanalysis(:);
BLoss_value_31to75 = BLoss_value_tmp(ssIters:ssIters2);
FrobLoss_value_tmp = Frobenius_norm_precision_finalanalysis(:);
FrobLoss_value_31to75 = FrobLoss_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2

ssIters = 76;
ssIters2 = 83;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 76;
ssIters2 = 83;

ELoss_value_tmp = entropy_loss_finalanalysis(:);
ELoss_value_76to83 = ELoss_value_tmp(ssIters:ssIters2);
BLoss_value_tmp = bounded_loss_finalanalysis(:);
BLoss_value_76to83 = BLoss_value_tmp(ssIters:ssIters2);
FrobLoss_value_tmp = Frobenius_norm_precision_finalanalysis(:);
FrobLoss_value_76to83 = FrobLoss_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2


ssIters = 84;
ssIters2 = 100;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 84;
ssIters2 = 100;

ELoss_value_tmp = entropy_loss_finalanalysis(:);
ELoss_value_84to100 = ELoss_value_tmp(ssIters:ssIters2);
BLoss_value_tmp = bounded_loss_finalanalysis(:);
BLoss_value_84to100 = BLoss_value_tmp(ssIters:ssIters2);
FrobLoss_value_tmp = Frobenius_norm_precision_finalanalysis(:);
FrobLoss_value_84to100 = FrobLoss_value_tmp(ssIters:ssIters2);

ELoss_value = [ELoss_value_1to25; ELoss_value_26to30; ELoss_value_31to75;...
    ELoss_value_76to83; ELoss_value_84to100];

BLoss_value = [BLoss_value_1to25; BLoss_value_26to30; BLoss_value_31to75;...
    BLoss_value_76to83; BLoss_value_84to100];

FrobLoss_value = [FrobLoss_value_1to25; FrobLoss_value_26to30; FrobLoss_value_31to75;...
    FrobLoss_value_76to83; FrobLoss_value_84to100];

ELoss_Type = repmat({'ELoss'}, [reps,1]);
BLoss_Type = repmat({'BLoss'}, [reps,1]);
FrobLoss_Type = repmat({'FrobLoss'}, [reps,1]);

Loss_value = [ELoss_value; BLoss_value; FrobLoss_value];
Loss_Type = [ELoss_Type; BLoss_Type; FrobLoss_Type];

[number_elements, ~] = size(Loss_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_twopercent_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_twopercent_Bsplines];


clearvars -except combine_tables


%AR4 p100 n500 Bsplines 

load('RankLikelihood_p100_n500_AR4_Bsplines_final.mat');

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

table_p100_n500_AR4_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR4_Bsplines];


clearvars -except combine_tables

%AR1 p100 n500 Bsplines 

load('RankLikelihood_p100_n500_AR1_Bsplines_final.mat');

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

table_p100_n500_AR1_Bsplines = table(Loss_value, Method, Loss_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR1_Bsplines];


%Try writing table as a csv file to read into R

writetable(combine_tables,'RankLikelihood_Boxplot_Loss_p100_n500.csv') 

