#
# Given a network and disease, find combinations of drugs that may be useful
#

library("igraph")
library("dplyr")

source("scripts/filterNodes.r")
source("scripts/filterEdges.r")
source("scripts/findDistance.r")

# Find pairs of drugs that may be useful
findDrugPairs <- function(graph, chosenDisease, order) {
    dist <- getDistFromDisease(graph, chosenDisease, order)

    drugNames <- row.names(dist[dist$diseaseDist < 0, ])
    effectiveness <- dist[dist$diseaseDist < 0, "diseaseDist", drop = FALSE]
    drugSep <- dist[drugNames, drugNames, drop = FALSE]

    separated <- which(drugSep > 0, arr.ind = TRUE)
    drugPairs <- data.frame(
        Drug1 = rownames(drugSep)[pmin(separated[, 1], separated[, 2])],
        Drug2 = colnames(drugSep)[pmax(separated[, 1], separated[, 2])],
        s_AB = unlist(drugSep[separated])
    ) %>% distinct(Drug1, Drug2, .keep_all = TRUE)

    return(list(drugPairs, effectiveness))
}

### Sanity Check
# findDrugPairs(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 1)




