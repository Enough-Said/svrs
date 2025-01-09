#
# A pipeline displaying how to use all the scripts and functions together to find
# drug combinations
#

source("analysisScripts/filterNodes.r")
source("analysisScripts/filterEdges.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")

################################################################################

chosenDisease <- "ALZHEIMER DISEASE 2"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.tissue(g)

drugDist <- findDrugDist(g, chosenDisease, 0, 0)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugDist, drugPairs, maxSize = 3)

print(drugCombs)

################################################################################

chosenDisease <- "Alzheimer's Disease"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.type(g)
newg <- v.filter.tissue(g, "hippocampal formation", minNTPM = 1)

drugDist <- findDrugDist(g, chosenDisease, 1, 1, overlapDist)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugDist, drugPairs, maxSize = 3)

print(drugCombs[[1]])

################################################################################

# Metformin in combination with temozolomide showed promising synergistic efficacy for treatment of glioblastoma
chosenDisease <- "Glioblastoma"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)

tissue <- fromJSON("clean/tissue.json")
tCutoff <- data.frame(
    tissue = unique(tissue$Tissue),
    min = rep(10, times = length(unique(tissue$Tissue))),
    max = rep(1000000, times = length(unique(tissue$Tissue)))
)
g <- e.filter.tissue(g, ntpmCutoff = tCutoff)

drugDist <- findDrugDist(g, chosenDisease, 1, 1)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugDist, drugPairs, maxSize = 2)

print(drugCombs)
