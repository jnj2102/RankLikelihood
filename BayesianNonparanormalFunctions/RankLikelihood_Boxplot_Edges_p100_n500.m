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
%Edges_Type = {'Specificity', 'Sensitivity', 'MCC'};
%Sparsity = {'AR4', 'Star', 'AR1', 'Percent'};

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_twopercent_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [table_p100_n500_twopercent_ranks];


clearvars -except combine_tables


%AR4 p100 n500 ranks 

load('RankLikelihood_p100_n500_AR4_ranks_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR4'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_AR4_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR4_ranks];


clearvars -except combine_tables

%AR1 p100 n500 ranks 

load('RankLikelihood_p100_n500_AR1_ranks_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR1'}, [number_elements,1]);
Method = repmat({'Ranks'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_AR1_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR1_ranks];


clearvars -except combine_tables

%twopercent p100 n500 Bsplines 
ssIters = 1;
ssIters2 = 25;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 1;
ssIters2 = 25;

SP_value_tmp = SP_matrix_finalanalysis(:);
SE_value_tmp = SE_matrix_finalanalysis(:);
MCC_value_tmp = MCC_matrix_finalanalysis(:);

SP_value_1to25 = SP_value_tmp(ssIters:ssIters2);
SE_value_1to25 = SE_value_tmp(ssIters:ssIters2);
MCC_value_1to25 = MCC_value_tmp(ssIters:ssIters2);
 
clear ssIters ssIters2

ssIters = 26;
ssIters2 = 30;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 26;
ssIters2 = 30;

SP_value_tmp = SP_matrix_finalanalysis(:);
SE_value_tmp = SE_matrix_finalanalysis(:);
MCC_value_tmp = MCC_matrix_finalanalysis(:);

SP_value_26to30 = SP_value_tmp(ssIters:ssIters2);
SE_value_26to30 = SE_value_tmp(ssIters:ssIters2);
MCC_value_26to30 = MCC_value_tmp(ssIters:ssIters2);


clear ssIters ssIters2

ssIters = 31;
ssIters2 = 75;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 31;
ssIters2 = 75;

SP_value_tmp = SP_matrix_finalanalysis(:);
SE_value_tmp = SE_matrix_finalanalysis(:);
MCC_value_tmp = MCC_matrix_finalanalysis(:);

SP_value_31to75 = SP_value_tmp(ssIters:ssIters2);
SE_value_31to75 = SE_value_tmp(ssIters:ssIters2);
MCC_value_31to75 = MCC_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2

ssIters = 76;
ssIters2 = 83;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 76;
ssIters2 = 83;

SP_value_tmp = SP_matrix_finalanalysis(:);
SE_value_tmp = SE_matrix_finalanalysis(:);
MCC_value_tmp = MCC_matrix_finalanalysis(:);

SP_value_76to83 = SP_value_tmp(ssIters:ssIters2);
SE_value_76to83 = SE_value_tmp(ssIters:ssIters2);
MCC_value_76to83 = MCC_value_tmp(ssIters:ssIters2);

clear ssIters ssIters2

ssIters = 84;
ssIters2 = 100;

%Load the first set
load(sprintf('RankLikelihood_p100_n500_twopercent_Bsplines_final_%dto%d.mat', ssIters, ssIters2));

ssIters = 84;
ssIters2 = 100;

SP_value_tmp = SP_matrix_finalanalysis(:);
SE_value_tmp = SE_matrix_finalanalysis(:);
MCC_value_tmp = MCC_matrix_finalanalysis(:);

SP_value_84to100 = SP_value_tmp(ssIters:ssIters2);
SE_value_84to100 = SE_value_tmp(ssIters:ssIters2);
MCC_value_84to100 = MCC_value_tmp(ssIters:ssIters2);

SP_value = [SP_value_1to25; SP_value_26to30; SP_value_31to75; SP_value_76to83;...
    SP_value_84to100];

SE_value = [SE_value_1to25; SE_value_26to30; SE_value_31to75; SE_value_76to83;...
    SE_value_84to100];

MCC_value = [MCC_value_1to25; MCC_value_26to30; MCC_value_31to75; MCC_value_76to83;...
    MCC_value_84to100];

SP_Type = repmat({'Specificity'}, [reps,1]);
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_twopercent_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_twopercent_Bsplines];


clearvars -except combine_tables


%AR4 p100 n500 Bsplines 

load('RankLikelihood_p100_n500_AR4_Bsplines_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR4'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_AR4_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR4_Bsplines];


clearvars -except combine_tables

%AR1 p100 n500 Bsplines 

load('RankLikelihood_p100_n500_AR1_Bsplines_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR1'}, [number_elements,1]);
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_AR1_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR1_Bsplines];



%Try writing table as a csv file to read into R

writetable(combine_tables, 'RankLikelihood_Boxplot_Edges_p100_n500.csv') 

