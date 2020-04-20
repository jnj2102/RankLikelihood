#Try doing a little graph analysis and compare the results between the graphs.  See if the analysis means something:
#So see if the node with max edges means something, the cluster of nodes mean something (community structure, clique).  
#See if they are similar or different between the graphs and try to relate to biology (what the nodes mean).
#Author: Jami Jackson Mulgrave
######################################################################################################################

install.packages("igraph")
install.packages("R.matlab")
installed.packages("BDgraph")
library(igraph)
library(R.matlab)
library(BDgraph)


#Need to read v7.3 mat files
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
library(rhdf5)

#read in the original data because we need the variable names


#Data from BDGraph package

data <- readMat('BDGraph_dataset.mat', package = "R.matlab") #Read in all the data

Realdata_variable_names = data$variableNamesGeneExpression
Realdata_variable_names_vector = matrix(unlist(Realdata_variable_names), ncol = 1)

####Read in the mat file for B-splines method####

Bsplines_adjmatrix <- h5read("C:/Users/jnjac/Documents/Ghoshal/Draft Papers/Paper_3/Real_Data/Paper_3/RealData_RankLikelihood_Bsplines.mat", "edge_matrix_finalanalysis")

#take the graph from the adjacency matrix

#not including the diagonal in the calculations 
Bsplines_graph = graph_from_adjacency_matrix(Bsplines_adjmatrix, mode = c("undirected"), diag = FALSE)

#find the degree

deg_Bsplines_graph <- degree(Bsplines_graph, mode="all")

#create a matrix of the degree and variable names

deg_Bsplines_graph_df = data.frame(deg_Bsplines_graph, Realdata_variable_names_vector)
colnames(deg_Bsplines_graph_df) <- c("Degree", "Variable")

#Order the degree

ordered_deg_Bsplines_graph_df = deg_Bsplines_graph_df[order(-deg_Bsplines_graph_df$Degree), ] #11, 10,10,9,7,7,6,6,6...

#each position corresponds to the node.

####Repeat for the Rank Likelihood method####

#Read in the mat file for B-splines method

Ranks_adjmatrix <- h5read("C:/Users/jnjac/Documents/Ghoshal/Draft Papers/Paper_3/Real_Data/Paper_3/RealData_RankLikelihood_ranks.mat", "edge_matrix_finalanalysis")

#take the graph from the adjacency matrix

#not including the diagonal in the calculations 
Ranks_graph = graph_from_adjacency_matrix(Ranks_adjmatrix, mode = c("undirected"), diag = FALSE)
#out of the top 6, 3 are the same for ranks and b-splines - take care with ties though.  No ties - list all top xx degrees
#these variables that match may mean something important - very connected.  
#Also compare this with the papers - See if the top nodes are in these table 2 and pg 27.

#find the degree

deg_Ranks_graph <- degree(Ranks_graph, mode="all")

#create a matrix of the degree and variable names

deg_Ranks_graph_df = data.frame(deg_Ranks_graph, Realdata_variable_names_vector)
colnames(deg_Ranks_graph_df) <- c("Degree", "Variable")

#Order the degree

ordered_deg_Ranks_graph_df = deg_Ranks_graph_df[order(-deg_Ranks_graph_df$Degree), ] 

#each position corresponds to the node.

####BDGraph Graph ####


#Read in the mat file for BDGraph method

BDGraph_adjmatrix <- readMat("C:/Users/jnjac/Documents/Ghoshal/Draft Papers/Paper_3/Real_Data/Paper_3/RealData_BDGraph.mat")

#convert to matrix

BDGraph_adjmatrix = matrix(unlist(BDGraph_adjmatrix), ncol = length(Realdata_variable_names_vector))
#take the graph from the adjacency matrix

#not including the diagonal in the calculations 
BDGraph_graph = graph_from_adjacency_matrix(BDGraph_adjmatrix, mode = c("undirected"), diag = FALSE)

#find the degree

deg_BDGraph_graph <- degree(BDGraph_graph, mode="all")

#create a matrix of the degree and variable names

deg_BDGraph_graph_df = data.frame(deg_BDGraph_graph, Realdata_variable_names_vector)
colnames(deg_BDGraph_graph_df) <- c("Degree", "Variable")

#Order the degree

ordered_deg_BDGraph_graph_df = deg_BDGraph_graph_df[order(-deg_BDGraph_graph_df$Degree), ] #shares 2 with bsplines

#each position corresponds to the node.


#Read in the mat file for Truncation method

Truncation_adjmatrix <- readMat("C:/Users/jnjac/Documents/Ghoshal/Draft Papers/Paper_3/Real_Data/Paper_3/RealData_truncation.mat")

Truncation_adjmatrix = matrix(unlist(Truncation_adjmatrix), ncol = length(Realdata_variable_names_vector))
#take the graph from the adjacency matrix

#not including the diagonal in the calculations 
Truncation_graph = graph_from_adjacency_matrix(Truncation_adjmatrix, mode = c("undirected"), diag = FALSE)

#find the degree

deg_Truncation_graph <- degree(Truncation_graph, mode="all")

#create a matrix of the degree and variable names

deg_Truncation_graph_df = data.frame(deg_Truncation_graph, Realdata_variable_names_vector)
colnames(deg_Truncation_graph_df) <- c("Degree", "Variable")

#Order the degree

ordered_deg_Truncation_graph_df = deg_Truncation_graph_df[order(-deg_Truncation_graph_df$Degree), ] #shares 2 with BDgraph

#each position corresponds to the node.



#find the intersection of these dataframes to see what they share.



#find community structure 


#find the clique
#cliques(Bsplines_graph)  #there are 240 - don't really want to go through these.

largest_cliques_Bsplines = largest_cliques(Bsplines_graph)  #cliques with max number of nodes  - this might be interesting.  
#Has 4 nodes in the cliques

largest_cliques_Ranks = largest_cliques(Ranks_graph) #has 5 nodes in the cliques

largest_cliques_BDGraph = largest_cliques(BDGraph_graph) #has 5 nodes in the cliques

largest_cliques_Truncation = largest_cliques(Truncation_graph) #has the most number of nodes in its cliques - 9

#not sure about cliques - nothing seems to match.

#Looking at the code - I don't think community detection will be useful either.


ceb_Bsplines <- cluster_edge_betweenness(Bsplines_graph) 
length(ceb_Bsplines)  #48 communities
membership(ceb_Bsplines)

#what nodes are part of the most frequent occurring community?
table(membership(ceb_Bsplines)) #community 28 comes up 16 times, community 15 comes up 13 ties.



ceb_Ranks <- cluster_edge_betweenness(Ranks_graph) 
length(ceb_Ranks)  #48 communities
membership(ceb_Ranks)

#what nodes are part of the most frequent occurring community?
table(membership(ceb_Ranks)) #community 26 comes up 10 times, community 29 comes up 9 times


######################################################################################################
#Transcripts from Joint High-Dimensional Bayesian Variable and Covariance Selection

#55 transcripts

transcripts_JointHDBVCS = c("GI_7019408-S", "GI_4504436-S", "GI_20302136-S", "GI_7661757-S", "GI_28610153-S",
                            "GI_4505888-A",  "GI_16554578-S", "GI_27754767-A","GI_41350202-S", "GI_20070269-S",
                            "GI_18379361-A", "GI_20373176-S", "GI_14211892-S", "GI_13514808-S", "GI_17981706-S",
                            "GI_33356559-S", "GI_33356162-S", "GI_40354211-S", "GI_27754767-I", "GI_9961355-S",
                            "GI_22027487-S", "hmm9615-S", "GI_13027804-S","GI_38569448-S", "GI_37537697-S",
                            "GI_34222299-S", "GI_28416938-S", "GI_18641371-S", "GI_21614524-S", "GI_37537705-I",
                            "GI_41190543-S", "GI_31652245-I", "GI_30795192-A", "GI_27482629-S", 
                            "GI_41197088-S", "GI_28557780-S", "GI_23510363-A",  "hmm3574-S", 
                            "GI_16159362-S", "GI_21389558-S", "GI_24308084-S", "GI_18641372-S", "GI_41190507-S",
                            "GI_37546969-S", "GI_5454143-S", "GI_27477086-S", "GI_18426974-S", "GI_4504410-S", 
                            "GI_27894333-S", "GI_19224662-S", "GI_19743804-A", "GI_4504700-S",
                            "GI_4502630-S")


length(transcripts_JointHDBVCS) #53

#there are only 53 - need to find 2 more in the graph but the picture is off because there are a bunch of Hs
#but the names have more than Hs.

transcripts_LGHD = c("GI_41197088-S", "GI_18641372-S", "GI_4504410-S", "GI_18426974-S", "GI_11095446-S", "GI_20302136-S
", "GI_7661757-S", "GI_41190507-S", "Hs.512137-S", "Hs.512124-S", "Hs.449605-S", "GI_37546969-S", "hmm3574-S",
                     "Hs.406489-S", "GI_27482629-S", "GI_13325059-S", "GI_19224662-S", "GI_23510363-A",
                     "GI_19743804-A", "GI_31377723-S", "GI_27894333-A", "GI_4502630-S", "GI_28557780-S",
                     "GI_27477086-S", "GI_27764881-S", "GI_5454143-S", "GI_16159362-S", "GI_24308084-S",
                     "GI_21389558-S", "GI_31652245-I" , "Hs.449609-S", "GI_41190543-S", "Hs.171273-S",
                     "GI_13027804-S", "hmm9615-S", "GI_37537697-S", "GI_34222299-S", "GI_22027487-S",
                     "GI_21614524-S", "GI_34915989-S", "GI_18379361-A", "GI_4504436-S", "GI_28610153-S",
                     "GI_20070269-S", "GI_14211892-S", "GI_33356559-S", "GI_33356162-S", "GI_20373176-S",
                     "GI_17981706-S", "GI_13514808-S", "GI_21464138-A", "GI_16554578-S", "GI_9961355-S",
                     "GI_41350202-S", "GI_38569448-S", "GI_27754767-A", "GI_27754767-I", "Hs.185140-S",
                     "GI_40354211-S", "GI_4505888-A")

length(transcripts_LGHD) #60 transcripts

duplicated(transcripts_LGHD) #no duplicates

#Bsplines

#compare the graphs themselves by only looking at nodes with at least 1 degree

ordered_deg_Bsplines_graph_df_degree1up = ordered_deg_Bsplines_graph_df[ordered_deg_Bsplines_graph_df$Degree >= 1, ]
dim(ordered_deg_Bsplines_graph_df_degree1up) #65 transcriipts

intersect_bsplines_JointHDBVCS_wholegraph =  intersect(ordered_deg_Bsplines_graph_df_degree1up$Variable, transcripts_JointHDBVCS)
length(intersect_bsplines_JointHDBVCS_wholegraph) #47 transcripts

#top 5 has ties.

#top degree is 11, then  two 10a

#compare the and top degree nodes  No truth.
# 
# ordered_deg_Bsplines_graph_df_degree5up = ordered_deg_Bsplines_graph_df[ordered_deg_Bsplines_graph_df$Degree >= 5, ]
# dim(ordered_deg_Bsplines_graph_df_degree5up) #12 transcripts
# 
# intersect_bsplines_JointHDBVCS_topnodes = intersect(ordered_deg_Bsplines_graph_df_degree5up$Variable, transcripts_JointHDBVCS)
# length(intersect_bsplines_JointHDBVCS_topnodes) #10 transcripts
# 
# #just looking at the top degree nodes is arbitrary.  I think just compare the edges - what do my graphs capture compared
#to these papers'?

#Ranks


#compare the graphs themselves by only looking at nodes with at least 1 degree

ordered_deg_Ranks_graph_df_degree1up = ordered_deg_Ranks_graph_df[ordered_deg_Ranks_graph_df$Degree >= 1, ]
dim(ordered_deg_Ranks_graph_df_degree1up) #92 transcriipts

intersect_Ranks_JointHDBVCS_wholegraph =  intersect(ordered_deg_Ranks_graph_df_degree1up$Variable, transcripts_JointHDBVCS)
length(intersect_Ranks_JointHDBVCS_wholegraph) #51 transcripts - so we captured a few more - almost all of them (plus more transcripts not there)

#top degree is 17, then two 16


#top 5 has ties.
#compare the and top degree nodes  No truth.

#ordered_deg_Ranks_graph_df_degree5up = ordered_deg_Ranks_graph_df[ordered_deg_Ranks_graph_df$Degree >= 5, ]
#dim(ordered_deg_Ranks_graph_df_degree5up) #44 transcripts

#intersect_Ranks_JointHDBVCS_topnodes = intersect(ordered_deg_Ranks_graph_df_degree5up$Variable, transcripts_JointHDBVCS)
#length(intersect_Ranks_JointHDBVCS_topnodes) #29 transcripts

#Truncation

#compare the graphs themselves by only looking at nodes with at least 1 degree

ordered_deg_Truncation_graph_df_degree1up = ordered_deg_Truncation_graph_df[ordered_deg_Truncation_graph_df$Degree >= 1, ]
dim(ordered_deg_Truncation_graph_df_degree1up) #96 transcriipts

intersect_Truncation_JointHDBVCS_wholegraph =  intersect(ordered_deg_Truncation_graph_df_degree1up$Variable, transcripts_JointHDBVCS)
length(intersect_Truncation_JointHDBVCS_wholegraph) # 52 transcripts

#could look at top 5 for truncation

#top 8?

#top degree is 21 (two of them).  Then three 19.


#compare the and top degree nodes  No truth.
# 
# ordered_deg_Truncation_graph_df_degree5up = ordered_deg_Truncation_graph_df[ordered_deg_Truncation_graph_df$Degree >= 5, ]
# dim(ordered_deg_Truncation_graph_df_degree5up) # 59 transcripts
# 
# intersect_Truncation_JointHDBVCS_topnodes = intersect(ordered_deg_Truncation_graph_df_degree5up$Variable, transcripts_JointHDBVCS)
# length(intersect_Truncation_JointHDBVCS_topnodes) # 40 transcripts
# 
# 
# #BDGraph

#compare the graphs themselves by only looking at nodes with at least 1 degree

ordered_deg_BDGraph_graph_df_degree1up = ordered_deg_BDGraph_graph_df[ordered_deg_BDGraph_graph_df$Degree >= 1, ]
dim(ordered_deg_BDGraph_graph_df_degree1up) #  100 transcriipts

intersect_BDGraph_JointHDBVCS_wholegraph =  intersect(ordered_deg_BDGraph_graph_df_degree1up$Variable, transcripts_JointHDBVCS)
length(intersect_BDGraph_JointHDBVCS_wholegraph) #  52 transcripts

#could look at top 5 for BDGraph

#top degree is 35, then 28, 27


#compare the and top degree nodes  No truth.
# 
# ordered_deg_BDGraph_graph_df_degree5up = ordered_deg_BDGraph_graph_df[ordered_deg_BDGraph_graph_df$Degree >= 5, ]
# dim(ordered_deg_BDGraph_graph_df_degree5up) #  100 transcripts
# 
# intersect_BDGraph_JointHDBVCS_topnodes = intersect(ordered_deg_BDGraph_graph_df_degree5up$Variable, transcripts_JointHDBVCS)
# length(intersect_BDGraph_JointHDBVCS_topnodes) #  52 transcripts