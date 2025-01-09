#
# Find topological distances between network nodes
#

library("igraph")
library("dplyr")
library("pbapply")
library("parallel")

# Find the average shortest distances between two sets of nodes
# Implements formula: s_AB = d_AB - (d_AA + d_BB)/2
# Where d_XY is the sum of shortest distances from nodes in X to any node in Y 
# and nodes in Y to any node in X, divided by the total number of nodes
topDist <- function(graph, nodeset1, nodeset2) {
    changeInfTo <- length(V(graph)) # May change if better method found
    
    d <- distances(graph, v = nodeset1, to = nodeset2, mode = "out")
    d[d == Inf] <- changeInfTo 
    d_AB <- sum(apply(d, 1, min)) + sum(apply(d, 2, min))
    d_AB <- d_AB / (nrow(d) + ncol(d))

    d <- distances(graph, v = nodeset1, to = nodeset1, mode = "out")
    diag(d) <- Inf # Remove the zeroes in diagonal
    d[d == Inf] <- changeInfTo
    d_AA <- sum(apply(d, 1, min))
    d_AA <- d_AA / nrow(d)

    d <- distances(graph, v = nodeset2, to = nodeset2, mode = "out")
    diag(d) <- Inf
    d[d == Inf] <- changeInfTo
    d_BB <- sum(apply(d, 1, min))
    d_BB <- d_BB / nrow(d)

    return(d_AB - (d_AA + d_BB)/2)
}

# Alternative distance function given by `shortest separation - overlap`
overlapDist <- function(graph, nodeset1, nodeset2) {
    changeInfTo <- length(V(graph))

    d <- min(distances(graph, v = nodeset1, to = nodeset2, mode = "out"))
    ans <- ifelse(d == Inf, changeInfTo, d - sum(nodeset1 %in% nodeset2))
    return(ans)
}


# Given a graph and disease name, finds nearby drugs and returns their topological distance values
# `minSep` and `maxSep` are the minimum and maximum separation, 
# representing the number of nodes that can appear between the disease and drugs
getDistFromDisease <- function(graph, d, minSep = 0, maxSep = 0, distFun = topDist) {
    v <- ego(graph, order = maxSep+2, nodes = d, mode = "all")[[1]]
    v <- v %m% ego(graph, order = minSep+1, nodes = d, mode = "all")[[1]]

    drugs <- v[v$type == "drug"]
    drug_names <- names(drugs)
    cat(paste("Found", length(drug_names), "drugs", "\n"))

    diseaseCluster <- neighbors(graph, d, mode = "out")
    drugCluster <- lapply(drugs, function(x) neighbors(graph, x, mode = "out"))

    cat("Calculating disease-drug distance\n")
    diseaseDist <- unlist(pblapply(drug_names,
        function(x) distFun(graph, diseaseCluster, drugCluster[[x]])))

    diseaseDist <- diseaseDist - maxSep
    drug_names <- drug_names[diseaseDist <= 0]
    diseaseDist <- diseaseDist[diseaseDist <= 0]

    cat(paste("Found", length(drug_names), "useful drugs\n"))
    if (length(drug_names) == 0) {
        return(data.frame())
    }

    cat("Calculating drug-drug distance\n")
    distM <- do.call(rbind, pblapply(drug_names, function(i) {
        mclapply(drug_names, 
            function(j) distFun(graph, drugCluster[[i]], drugCluster[[j]])
        )
    }))
    dimnames(distM) <- list(drug_names, drug_names)

    return(data.frame(cbind(diseaseDist, distM)))
}

### Sanity Check
g <- readRDS("clean/baseGraph.rds")
getDistFromDisease(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 0, 0)
getDistFromDisease(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 0, 0, overlapDist)
