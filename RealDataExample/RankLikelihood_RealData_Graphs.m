%Code to create the estimated graphs for the Real Data application
%
%Author: Jami Jackson Mulgrave
clear;

addpath('C:/Users/jnjac/Documents/MATLAB/paul-kassebaum-mathworks-circularGraph-3a7926b');

load('RealData_truncation.mat');
load('RealData_BDGraph.mat');

%ranks
load('RealData_RankLikelihood_ranks.mat');

myColorMap = repmat([0 0 0], [p,1]);



%Bayesian nonparanormal graph
%edge_matrix_finalanalysis(logical(eye(size(edge_matrix_finalanalysis)))) = 0;

edge_matrix_finalanalysis_ranks = double(edge_matrix_finalanalysis);

%call the circular graph

figure;

H = circularGraph(edge_matrix_finalanalysis_ranks, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes_ranks.fig')

figs = openfig('RealData_circularGraph_Bayes_ranks.fig');
   saveas(figs, 'RealData_circularGraph_Bayes_ranks.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere


clear edge_matrix_finalanalysis

%read in the Bsplines

%Bsplines
load('RealData_RankLikelihood_Bsplines.mat');


%Bayesian nonparanormal graph
%edge_matrix_finalanalysis(logical(eye(size(edge_matrix_finalanalysis)))) = 0;

edge_matrix_finalanalysis_Bsplines = double(edge_matrix_finalanalysis);

%call the circular graph

figure;

H = circularGraph(edge_matrix_finalanalysis_Bsplines, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes_Bsplines.fig')

figs = openfig('RealData_circularGraph_Bayes_Bsplines.fig');
   saveas(figs, 'RealData_circularGraph_Bayes_Bsplines.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere




%Read in the frequentist truncation model


figure;

truncation_graph = circularGraph(edgeMat_glasso, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_truncation.fig')

figs = openfig('RealData_circularGraph_Frequentist_truncation.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_truncation.png');
   

   %BDgraph 
   
   figure;

BDGraph_graph = circularGraph(edgeMat_BDGraph, 'Colormap',myColorMap);

savefig('RealData_circularGraph_BDGraph.fig')

figs = openfig('RealData_circularGraph_BDGraph.fig');
   saveas(figs, 'RealData_circularGraph_BDGraph.png');

%Read in the frequentist EBIC nonpararnormal model

   
   %How many different edges per graph?
   

indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 
  sum_edges_bayes_ranks =  sum(edge_matrix_finalanalysis_ranks(upperind) == 1);

    sum_edges_bayes_Bsplines =  sum(edge_matrix_finalanalysis_Bsplines(upperind) == 1);

    sum_edges_truncation =  sum(edgeMat_glasso(upperind) == 1);

                    sum_edges_BDGraph =  sum(edgeMat_BDGraph(upperind) == 1);
