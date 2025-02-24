#
# Get useful statistics for reporting
# 

library("igraph")

g <- readRDS("clean/baseGraph.rds")
length(V(g)[[type == "drug"]])
length(V(g)[[type == "human-protein"]])
length(V(g)[[type == "disease"]])
