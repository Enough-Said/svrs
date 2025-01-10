#
# Create a visualisation of the subgraph formed by the disease and srug combination
#

set.seed(123)
library("igraph")

plotCombination <- function(graph, disease, drugs, additionalNodes = NULL, diseaseOrder = 1, otherOrder = 1) {
    type_colors <- c("human-protein" = "lightblue", 
        "drug" = "green", 
        "disease" = "red")
    
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

    plot(
        as_undirected(subgraph),
        vertex.color = type_colors[V(subgraph)$type],
        main = paste("Subgraph of", disease, " and Given Drugs")
    )
}

getNearbyNodes <- function(graph, order, node) {
    return(names(unlist(ego(graph, order = order, nodes = node, mode = "out"))))
}


### Sanity Check
# chosenDisease <- "ALZHEIMER DISEASE 2"
# plotCombination(g, chosenDisease, unlist(drugCombs[[2]][[2]]),
#     c(
#         getNearbyNodes(g, 1, "CHEMBL463981"),
#         getNearbyNodes(g, 1, "CHEMBL3545261")
#     )
# )
