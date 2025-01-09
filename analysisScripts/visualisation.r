#
# Create a visualisation of the subgraph formed by the disease and srug combination
#

library("igraph")

plotCombination <- function(graph, disease, drugs, additionalNodes = NULL) {
    type_colors <- c("human-protein" = "lightblue", 
        "drug" = "green", 
        "disease" = "red")
    
    subgraph <- induced_subgraph(graph, unlist(ego(
        graph, 
        order = 1, 
        nodes = c(disease, drugs, additionalNodes),
        mode = "out"
    )))

    plot(
        subgraph,
        vertex.color = type_colors[V(subgraph)$type],
        main = "Subgraph with Vertices Colored by 'Type'"
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
