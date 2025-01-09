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

drugDist <- getDistFromDisease(g, chosenDisease, 0, 1)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugPairs, maxSize = 2)

print(drugCombs[[2]])

################################################################################

chosenDisease <- "Alzheimer's Disease"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.type(g)
newg <- v.filter.tissue(g, "hippocampal formation", minNTPM = 1)

drugDist <- getDistFromDisease(g, chosenDisease, 1, 1)
drugPairs <- findDrugPairs(g, drugDist)
drugCombs <- findDrugCombinations(g, drugPairs, maxSize = 3)

print(drugCombs[[1]])

################################################################################
