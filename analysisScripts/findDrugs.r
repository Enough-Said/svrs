#
# Given a network and disease, find combinations of drugs that may be useful
#

library("igraph")
library("dplyr")

source("analysisScripts/findDistance.r")

# Find all pairs of drugs that may be useful for the given disease
# Must provide a graph and all distances for each drug from the disease and from each other
# `minSep` is an optional parameter, and only drug pairs with separation strictly greater than
# it will be included
findDrugPairs <- function(graph, dist, minSep = 0) {
    drugNames <- row.names(dist[dist$diseaseDist <= 0, ])
    drugSep <- dist[drugNames, drugNames, drop = FALSE]

    separated <- which(drugSep > minSep, arr.ind = TRUE)
    drugPairs <- data.frame(
        Drug1 = rownames(drugSep)[pmin(separated[, 1], separated[, 2])],
        Drug2 = colnames(drugSep)[pmax(separated[, 1], separated[, 2])],
        dist = unlist(drugSep[separated])
    ) %>% distinct(Drug1, Drug2, .keep_all = TRUE)

    return(drugPairs)
}

# Find all combinations of drugs that can be useful in the chosen disease
# Can give a precomputed distance or precomputed drug pairs list
findDrugCombinations <- function(graph, drugDist, drugPairs, maxSize = -1) {
    items <- row.names(drugDist)

    drugPairs <- apply(drugPairs[, 1:2], 1, function(x) list(x[[1]], x[[2]]))
    pairSet <- unique(unlist(lapply(drugPairs, paste, collapse = "_")))
    isValidPair <- function(x, y) {
        paste(c(x, y), collapse = "_") %in% pairSet ||
        paste(c(y, x), collapse = "_") %in% pairSet
    }

    dp <- list()
    dp[[1]] <- lapply(items, function(item) c(item))

    maxSize <- ifelse(maxSize == -1, length(items), maxSize)
    for (size in 2:maxSize) {
        dp[[size]] <- list()

        for (comb in dp[[size - 1]]) {
            for (item in items) {
                if (item %in% comb) next

                # Check validity: item must form valid pairs with all items in the combination
                if (all(sapply(comb, isValidPair, y = item))) {
                    new_comb <- sort(c(comb, item))
                    
                    if (!any(sapply(dp[[size]], function(x) identical(x, new_comb)))) {
                        dp[[size]] <- c(dp[[size]], list(new_comb))
                    }
                }
            }
        }
    }

    return(dp)
}

### Sanity Check
# dist <- getDistFromDisease(g, "ALZHEIMER DISEASE 2", 0, 1)
# pairsAD2 <- findDrugPairs(g, dist)
# combs <- findDrugCombinations(g, dist, pairsAD2)

# print(combs)
