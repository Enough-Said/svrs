#
# Takes the base graph and filter nodes based on tissue and cell-line data.
# This file contains a series of functions for filtering the graph.
# 

library("jsonlite")
library("igraph")

g <- readRDS("clean/baseGraph.rds")
tissue <- fromJSON("clean/tissue.json")
cellLine <- read.delim("clean/cellLine.tsv")
cellType <- fromJSON("clean/cellType.json")
subcell <- fromJSON("clean/subcell.json")

# Given a dataframe containing the columns "Gene" and "nTPM", return a list of 
# genes where the nTPM is within the range minNTPM - maxNTPM inclusive. 
# vInG must be a list of vertices in the graph, found through `V(graph)`.
#
# Note that this returns only human PPI that are within the range given.
tpmFilter <- function(df, minNTPM, maxNTPM, vInG) {
    df <- df[df$nTPM >= minNTPM, ]
    df <- df[df$nTPM <= maxNTPM, ]

    df <- df[df$Gene %in% vInG$name, ]
    return(which(vInG$name %in% df$Gene))
}

# Filter the given graph based on the tissue t.
# Returns a subgraph containing only proteins with nTPM between the given values 
# for the tissue, as well as non-human-protein nodes. 
#
# nTPM = Normalised Transcripts per Million
v.filter.tissue <- function(graph, t, minNTPM = 1, maxNTPM = 1000000) {
    filteredGenes <- tpmFilter(tissue[tissue$Tissue == t, ], minNTPM, maxNTPM, 
                               V(graph))
    notPPI <- V(graph)$name[V(graph)$type != "human-protein"]
    notPPI <- which(V(graph)$name %in% notPPI)

    subgraph <- induced_subgraph(graph, vids = c(filteredGenes, notPPI))
    subgraph <- delete_vertices(subgraph, V(subgraph)[degree(subgraph) == 0])
    return(subgraph)
}

# Filter the given graph based on the cell type t.
# Returns a subgraph containing only proteins with nTPM between the given values 
# for the cell type, as well as non-human-protein nodes. 
#
# nTPM = Normalised Transcripts per Million
v.filter.type <- function(graph, t, minNTPM = 1, maxNTPM = 1000000) {
    filteredGenes <- tpmFilter(cellType[cellType$Cell.type == t, ], minNTPM, 
                               maxNTPM, V(graph))
    notPPI <- V(graph)$name[V(graph)$type != "human-protein"]
    notPPI <- which(V(graph)$name %in% notPPI)

    subgraph <- induced_subgraph(graph, vids = c(filteredGenes, notPPI))
    subgraph <- delete_vertices(subgraph, V(subgraph)[degree(subgraph) == 0])
    return(subgraph)
}

# Filter the given graph based on the cell line `cl`.
# Returns a subgraph containing only proteins with nTPM between the given values 
# for the cell line, as well as non-human-protein nodes. 
#
# nTPM = Normalised Transcripts per Million
v.filter.line <- function(graph, cl, minNTPM = 1, maxNTPM = 1000000) {
    filteredGenes <- tpmFilter(cellLine[cellLine$Cell.line == cl, ], minNTPM, 
                               maxNTPM, V(graph))
    notPPI <- V(graph)$name[V(graph)$type != "human-protein"]
    notPPI <- which(V(graph)$name %in% notPPI)

    subgraph <- induced_subgraph(graph, vids = c(filteredGenes, notPPI))
    subgraph <- delete_vertices(subgraph, V(subgraph)[degree(subgraph) == 0])
    return(subgraph)
}

# Filter the given graph based on subcellular location `l`. Returns a subgraph
# containing only genes expressed in the given location as well as 
# non-human-protein nodes.
v.filter.subcell <- function(graph, l) {
    filteredGenes <- subcell[sapply(
        subcell$All.location, function(x) l %in% x), ]$Gene

    filteredGenes <- filteredGenes[filteredGenes %in% V(graph)$name]
    filteredGenes <- which(V(graph)$name %in% filteredGenes)

    notPPI <- V(graph)$name[V(graph)$type != "human-protein"]
    notPPI <- which(V(graph)$name %in% notPPI)

    subgraph <- induced_subgraph(graph, vids = c(filteredGenes, notPPI))
    subgraph <- delete_vertices(subgraph, V(subgraph)[degree(subgraph) == 0])
    return(subgraph)
}




### Sanity Checks
# sg <- v.filter.tissue(g, "stomach", 0, 50)
# head(V(sg))
# 
# sg <- v.filter.type(g, "Adipocytes", 10, 20)
# head(V(sg))
# 
# head(cellLine[cellLine$Cell.line == "143B", ], n = 30)
# sg <- v.filter.line(g, "143B", 0, 50)
# head(V(sg))
# 
# sg <- v.filter.subcell(g, "Cytosol")
# head(V(sg))

