#
# Create a visualisation of the subgraph formed by the disease and srug combination
#

set.seed(123)
library("igraph")
library("networkD3")

plotCombination <- function(graph, disease, drugs, additionalNodes = NULL, 
    diseaseOrder = 1, otherOrder = 1, interactive = TRUE) {
    diseaseNodes <- unlist(ego(
        graph, 
        order = diseaseOrder, 
        nodes = disease,
        mode = "out"
    ))

    otherNodes <- unlist(ego(
        graph, 
        order = otherOrder, 
        nodes = c(drugs, additionalNodes),
        mode = "out"
    ))

    subgraph <- induced_subgraph(graph, c(diseaseNodes, otherNodes))

    if (!interactive) {
        type_colors <- c("human-protein" = "lightblue", 
            "drug" = "green", 
            "disease" = "red")
        plot(
            as_undirected(subgraph),
            vertex.color = type_colors[V(subgraph)$type],
            main = paste("Subgraph of", disease, " and Given Drugs")
        )
        return()
    }

    subgraphND3 <- igraph_to_networkD3(subgraph, group = V(subgraph)$type)
    D3plot <- forceNetwork(Links = subgraphND3$links, Nodes = subgraphND3$nodes, 
        Source = "source", Target = "target", 
        NodeID = "name", Group = "group", fontSize = 20, opacity = 1, 
        zoom = TRUE, colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
    
    print(D3plot)
}

getNearbyNodes <- function(graph, order, node) {
    return(names(unlist(ego(graph, order = order, nodes = node, mode = "out"))))
}


### Sanity Check
chosenDisease <- "ALZHEIMER DISEASE 2"
g <- readRDS("clean/baseGraph.rds")
drugDist <- findDrugDist(g, chosenDisease, 1, 1, overlapDist)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugDist, drugPairs, maxSize = 3)

plotCombination(g, chosenDisease, unlist(drugCombs[[3]][[1]]))
