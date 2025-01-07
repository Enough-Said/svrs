#
# Given a network and disease, find combinations of drugs that may be useful
#

library("igraph")
library("dplyr")

source("analysisScripts/findDistance.r")

# Find all pairs of drugs that may be useful for the given disease
# You can precompute distance for significant speed up if this function is repeated
# This also enables the use of other distance measures
findDrugPairs <- function(graph, chosenDisease, order, precompDist) {
    if (missing(precompDist)) {
        dist <- getDistFromDisease(graph, chosenDisease, order)
        if (nrow(dist) == 0) {
            message("No useful drugs found")
            return()
        }
    } else {
        dist <- precompDist
    }

    drugNames <- row.names(dist[dist$diseaseDist < 0, ])
    effectiveness <- dist[dist$diseaseDist < 0, "diseaseDist", drop = FALSE]
    drugSep <- dist[drugNames, drugNames, drop = FALSE]

    separated <- which(drugSep >= 0, arr.ind = TRUE)
    drugPairs <- data.frame(
        Drug1 = rownames(drugSep)[pmin(separated[, 1], separated[, 2])],
        Drug2 = colnames(drugSep)[pmax(separated[, 1], separated[, 2])],
        s_AB = unlist(drugSep[separated])
    ) %>% distinct(Drug1, Drug2, .keep_all = TRUE)

    return(list(drugPairs, effectiveness))
}

# Find all combinations of drugs that can be useful in the chosen disease
# Can give a precomputed distance or precomputed drug pairs list
findDrugCombinations <- function(graph, chosenDisease, order, precompDist, precompPairs, maxSize = -1) {
    dist <- NULL
    drugPairs <- NULL

    if (missing(precompPairs)) {
        if (missing(precompDist)) {
            drugPairs <- findDrugPairs(graph, chosenDisease, order)
        } else {
            dist <- precompDist
            drugPairs <- findDrugPairs(graph, chosenDisease, order, dist)
        }
    } else {
        drugPairs <- precompPairs
    }

    effectiveness <- drugPairs[[2]]
    drugPairs <- apply(drugPairs[[1]][, 1:2], 1, function(x) list(x[[1]], x[[2]]))
    pairSet <- unique(unlist(lapply(drugPairs, paste, collapse = "_")))
    isValidPair <- function(x, y) {
        paste(c(x, y), collapse = "_") %in% pairSet ||
        paste(c(y, x), collapse = "_") %in% pairSet
    }

    dp <- list()
    dp[[1]] <- lapply(row.names(effectiveness), function(item) c(item))

    maxSize <- ifelse(maxSize == -1, length(row.names(effectiveness)), maxSize)
    items <- row.names(effectiveness)
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

    return(list(dp, effectiveness))
}

### Sanity Check
# dist <- getDistFromDisease(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 1)
# pairsAMD1 <- findDrugPairs(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 1, dist)
# combs <- findDrugCombinations(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 1, precompPairs = pairsAMD1)

# print(combs)
