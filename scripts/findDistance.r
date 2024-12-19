#
# Find topological distances between network nodes
#

library("igraph")
library("dplyr")
library("pbapply")
library("parallel")

# Find the average shortest distances between two sets of nodes
# Implements formula: s_AB = d_AB - (d_AA + d_BB)/2
# Where d_XY is the sum of shortest distances from nodes in X to any node in Y and nodes in Y to any node in X, divided by the total number of nodes
top.dist <- function(graph, nodeset1, nodeset2) {
    changeInfTo <- length(V(graph)) # May change if better method found
    
    d <- distances(graph, v = nodeset1, to = nodeset2, mode = "out")
    d[d == Inf] <- changeInfTo 
    d_AB <- sum(apply(d, 1, min)) + sum(apply(d, 2, min))
    d_AB <- d_AB / (nrow(d) + ncol(d))

    d <- distances(graph, v = nodeset1, to = nodeset1, mode = "out")
    diag(d) <- max(d) # Remove the zeroes in diagonal
    d[d == Inf] <- changeInfTo
    d_AA <- sum(apply(d, 1, min))
    d_AA <- d_AA / nrow(d)

    d <- distances(graph, v = nodeset2, to = nodeset2, mode = "out")
    diag(d) <- max(d)
    d[d == Inf] <- changeInfTo
    d_BB <- sum(apply(d, 1, min))
    d_BB <- d_BB / nrow(d)

    return(d_AB - (d_AA + d_BB)/2)
}

# Given a graph and disease name, returns topological distance values
# `order` is the maximum amount of nodes that can appear between disease and drug proteins
getDistFromDisease <- function(graph, d, order) {
    v <- ego(graph, order = order+2, nodes = d, mode = "all")[[1]]
    drugs <- v[v$type == "drug"]
    drug_names <- names(drugs)

    diseaseCluster <- neighbors(graph, d, mode = "out")
    drugCluster <- lapply(drugs, function(x) neighbors(graph, x, mode = "out"))

    dist <- data.frame(
        diseaseDist = unlist(lapply(drug_names,
            function(x) top.dist(graph, diseaseCluster, drugCluster[[x]])
        ))
    )

    distM <- do.call(rbind, pblapply(seq_along(drug_names), function(i) {
        mclapply(seq_along(drug_names), function(j) {
            top.dist(graph, drugCluster[[i]], drugCluster[[j]])
        })
    }))
    dimnames(distM) <- list(drug_names, drug_names)

    return(cbind(dist, distM))
}

### Sanity Check
# getDistFromDisease(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 0)
