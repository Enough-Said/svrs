#
# Testing the effectiveness of the network-based repositioning method using known
# and verified pairwise drug combinations
#

library("jsonlite")
library("igraph")
library("tidyr")
library("dplyr")
library("Metrics")
library("pbapply")

source("analysisScripts/filterNodes.r")
source("analysisScripts/filterEdges.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")

source("analysisScripts/visualisation.r")

knownCombs <- fromJSON("clean/verifiedCombinations.json")

checkDrugDist <- function(graph, drugPairs, distFunc) {
    unlist(apply(drugPairs, MARGIN = 1, function(comb) {
        nodeset1 <- neighbors(graph, comb[[1]], mode = "out")
        nodeset2 <- neighbors(graph, comb[[2]], mode = "out")
        return(distFunc(graph, nodeset1, nodeset2))
    }))
}

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

##### Find predicted drug pairs that display complementary exposure
knownCombs$topDist <- checkDrugDist(g, knownCombs, topDist)
knownCombs$overlapDist <- checkDrugDist(g, knownCombs, overlapDist)

nrow(knownCombs) # 595 drug combinations
knownCombs %>% filter(overlapDist >= 0) %>% nrow() # 505 drug combinations
knownCombs %>% filter(topDist >= 0) %>% nrow() # 352 drug combinations


##### We can now repeat but with filtered edges
gFil <- readRDS("clean/baseGraph.rds")
gFil <- e.filter.subcell(gFil)
gFil <- e.filter.tissue(gFil)
gFil <- e.filter.type(gFil)

filtered <- knownCombs
filtered$topDist <- checkDrugDist(gFil, filtered, topDist)
filtered$overlapDist <- checkDrugDist(gFil, filtered, overlapDist)

nrow(filtered) # 595 drug combinations
filtered %>% filter(overlapDist >= 0) %>% nrow() # 505 drug combinations
filtered %>% filter(topDist >= 0) %>% nrow() # 462 drug combinations

##### We can check for false positives using randomm drug pairs
# Utilising method from https://www.nature.com/articles/s41467-019-09186-x

findAUC <- function(graph, drugPairs) {
    genPairs <- t(replicate(nrow(drugPairs), 
        sample(V(graph)[V(graph)$type == "drug"]$name, size = 2, replace = FALSE)))
    
    genPairsTopDist <- checkDrugDist(graph, drugPairs, topDist)

    finalAUC <- auc(c(rep(TRUE, nrow(drugPairs)), rep(FALSE, nrow(drugPairs))), 
        c(drugPairs$topDist >= 0, genPairsTopDist < 0))
    
    return(finalAUC)
}

noFiltMean <- mean(pbreplicate(100, findAUC(g, knownCombs))) # 0.5915966
filtMean <- mean(pbreplicate(100, findAUC(gFil, filtered)))   # 
