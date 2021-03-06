%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the RankLikelihood paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%tenPercent p25 n50 ranks 

load('RankLikelihood_p25_n50_tenpercent_ranks_final.mat');

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

table_p25_n50_tenpercent_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [table_p25_n50_tenpercent_ranks];


clearvars -except combine_tables


%AR4 p25 n50 ranks 

load('RankLikelihood_p25_n50_AR4_ranks_final.mat');

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

table_p25_n50_AR4_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n50_AR4_ranks];


clearvars -except combine_tables

%AR1 p25 n50 ranks 

load('RankLikelihood_p25_n50_AR1_ranks_final.mat');

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

table_p25_n50_AR1_ranks = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n50_AR1_ranks];


clearvars -except combine_tables

%tenPercent p25 n50 Bsplines 

load('RankLikelihood_p25_n50_tenpercent_Bsplines_final.mat');


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
Method = repmat({'Bsplines'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p25_n50_tenpercent_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n50_tenpercent_Bsplines];


clearvars -except combine_tables


%AR4 p25 n50 Bsplines 

load('RankLikelihood_p25_n50_AR4_Bsplines_final.mat');

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

table_p25_n50_AR4_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n50_AR4_Bsplines];


clearvars -except combine_tables

%AR1 p25 n50 Bsplines 

load('RankLikelihood_p25_n50_AR1_Bsplines_final.mat');

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

table_p25_n50_AR1_Bsplines = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p25_n50_AR1_Bsplines];



%Try writing table as a csv file to read into R

writetable(combine_tables,'RankLikelihood_Boxplot_Edges_p25_n50.csv') 

