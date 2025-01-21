#
# Testing the effectiveness of the network-based repositioning method using known
# and verified pairwise drug combinations
#

library("jsonlite")
library("igraph")
library("tidyr")
library("dplyr")

source("analysisScripts/filterNodes.r")
source("analysisScripts/filterEdges.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")

source("analysisScripts/visualisation.r")

knownCombs <- fromJSON("clean/verifiedCombinations.json")

# Since the drug combinations do not say what disease they are for, 
# we will find the closest diseases
g <- readRDS("clean/baseGraph.rds")
knownCombs$diseases <- apply(knownCombs, MARGIN = 1, function(drugs) {
    if (!(drugs[[1]] %in% V(g)$name && drugs[[2]] %in% V(g)$name)) {
        return(NULL)
    }
    dis <- ego(g, order = 2, nodes = c(drugs[[1]], drugs[[2]]), mode = "all")
    dis <- dis[[1]][dis[[1]] %in% dis[[2]]]
    dis <- dis[dis$type == "disease"]
    return(dis$name)
})

lengths <- sapply(knownCombs$diseases, function(x) length(x))
knownCombs <- knownCombs[lengths != 0, ]

# We can check the distances between drugs to hunt for complimentary exposure
knownCombs$topDist <- unlist(apply(knownCombs, MARGIN = 1, function(comb) {
    nodeset1 <- neighbors(g, comb[[1]], mode = "out")
    nodeset2 <- neighbors(g, comb[[2]], mode = "out")
    
    return(topDist(g, nodeset1, nodeset2))
}))

knownCombs$overlapDist <- unlist(apply(knownCombs, MARGIN = 1, function(comb) {
    nodeset1 <- neighbors(g, comb[[1]], mode = "out")
    nodeset2 <- neighbors(g, comb[[2]], mode = "out")
    
    return(overlapDist(g, nodeset1, nodeset2))
}))

nrow(knownCombs) # 595 drug combinations
knownCombs %>% filter(overlapDist >= 0) %>% nrow() # 505 drug combinations
knownCombs %>% filter(topDist >= 0) %>% nrow() # 352 drug combinations


### We can now repeat but with filtered edges
gFil <- readRDS("clean/baseGraph.rds")
gFil <- e.filter.subcell(gFil)
gFil <- e.filter.tissue(gFil)
gFil <- e.filter.type(gFil)

filtered <- knownCombs
filtered$topDist <- unlist(apply(filtered, MARGIN = 1, function(comb) {
    nodeset1 <- neighbors(gFil, comb[[1]], mode = "out")
    nodeset2 <- neighbors(gFil, comb[[2]], mode = "out")
    
    return(topDist(gFil, nodeset1, nodeset2))
}))

filtered$overlapDist <- unlist(apply(filtered, MARGIN = 1, function(comb) {
    nodeset1 <- neighbors(gFil, comb[[1]], mode = "out")
    nodeset2 <- neighbors(gFil, comb[[2]], mode = "out")
    
    return(overlapDist(gFil, nodeset1, nodeset2))
}))

nrow(filtered) # 595 drug combinations
filtered %>% filter(overlapDist >= 0) %>% nrow() # 505 drug combinations
filtered %>% filter(topDist >= 0) %>% nrow() # 462 drug combinations
