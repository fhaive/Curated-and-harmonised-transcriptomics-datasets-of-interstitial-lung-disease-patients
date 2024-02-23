rm(list=ls())

setwd("/path")
source("/path/INfORM_functions.R")

args = commandArgs(trailingOnly=TRUE)

file_path=args[1]


generatematrices=get_ranked_consensus_matrix(gx_table = read.table(file_path, sep="\t"), iMethods = c("clr"),
                                             iEst = c("pearson"),
                                             iDisc=c("none"), ncores = 30, debug_output = TRUE, updateProgress = TRUE)


#Parse ranked matrix and get bin_mat and edge_rank
# Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
# parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.

rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                      mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)

conGraph <- get_iGraph(rankMat.parsed$bin_mat)
saveRDS(conGraph, file="network.rds")


