#
# Takes the base graph, and removes impossible interactions
#

# We use variables and functions defined in this script
source("scripts/filterNodes.r")
library("dplyr")

### Credit to Jonno Bourne for union function
# The union of two or more graphs are created. 
# The graphs may have identical or overlapping vertex sets.
union2 <- function(g1, g2) {
    #Internal function that cleans the names of a given attribute
    CleanNames <- function(g, target) {
        # Get target names, find names with "_x" and remove
        gNames <- parse(text = (paste0(target,"_attr_names(g)"))) %>% eval 
        AttrNeedsCleaning <- grepl("(_\\d)$", gNames )
        StemName <- gsub("(_\\d)$", "", gNames)

        #replace attribute name for all attributes
        NewnNames <- unique(StemName[AttrNeedsCleaning])
        for (i in NewnNames) {
            attr1 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
            attr2 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval

            g <- parse(text = (paste0("set_",target,"_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))"))) %>%
                        eval

            g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
            g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval
        }
        return(g)
    }

    g <- igraph::union(g1, g2) 
    for (i in c("graph", "edge", "vertex")) {
        g <- CleanNames(g, i)
    }

    return(g)
}


# Remove human PPI that are physically separated by subcellular location
e.filter.subcell <- function(graph) {
    # Create each subgraph and combine them
    subcell.locations <- unique(unlist(subcell$All.location))
    newg <- Reduce(union2, 
        lapply(subcell.locations, function(l) {
            print(paste("Done", x[1]))
            v.filter.subcell(graph, l)
        }))
    return(newg)
}

# Only keep human PPI if both proteins are expressed in a tissue above a threshold
# Requires an igraph object and a dataframe with three columns: 
# Tissue, minimum nTPM and maximum nTPM
e.filter.tissue <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            tissue = unique(tissue$Tissue),
            min = rep(0.1, times = length(unique(tissue$Tissue))),
            max = rep(1000000, times = length(unique(tissue$Tissue)))
        )
    }

    newg <- Reduce(union2, 
        apply(ntpmCutoff, 1, function(x) {
            print(paste("Done", x[1]))
            v.filter.tissue(graph, x[1], x[2], x[3])
        }))

    return(newg)
}

# Only keep human PPI if both proteins are expressed in a cell type above a threshold
# Requires an igraph object and a dataframe with three columns: 
# Cell type, minimum nTPM and maximum nTPM
e.filter.type <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            cellType = unique(cellType$Cell.type),
            min = rep(0.1, times = length(unique(cellType$Cell.type))),
            max = rep(1000000, times = length(unique(cellType$Cell.type)))
        )
    }

    newg <- Reduce(union2, 
        apply(ntpmCutoff, 1, function(x) {
            print(paste("Done", x[1]))
            v.filter.type(graph, x[1], x[2], x[3])
        }))

    return(newg)
}

# Only keep human PPI if both proteins are expressed in a cell line above a threshold
# Requires an igraph object and a dataframe with three columns: 
# Cell line, minimum nTPM and maximum nTPM
e.filter.line <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            line = unique(cellLine$Cell.line),
            min = rep(0.1, times = length(unique(cellLine$Cell.line))),
            max = rep(1000000, times = length(unique(cellLine$Cell.line)))
        )
    }

    newg <- Reduce(union2, 
        apply(ntpmCutoff, 1, function(x) {
            print(paste("Done", x[1]))
            v.filter.line(graph, x[1], x[2], x[3])
        }))

    return(newg)
}

### Sanity Checks
# head(V(g)) # 32946
# head(E(g)) # 645539

# newg <- e.filter.subcell(g)
# head(V(newg)) # 24271
# head(E(newg)) # 440656

# # > 1200 cell lines; Very slow
# newg <- e.filter.line(g)
# head(V(newg))
# head(E(newg))

# newg <- e.filter.type(g)
# head(V(newg)) # 29790
# head(E(newg)) # 602840

# newg <- e.filter.tissue(g)
# head(V(newg)) # 29606
# head(E(newg)) # 596438
