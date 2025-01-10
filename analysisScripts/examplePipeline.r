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

# Verifying a known drug combination
g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.tissue(g)

chosenDisease <- "Non-small cell lung cancer stage I"
combo <- c("CHEMBL217092", "CHEMBL601719")

# Check the graph has all drugs
V(g)[combo]

checkDrugComb(g, chosenDisease, combo)
checkDrugComb(g, chosenDisease, combo, overlapDist)
plotCombination(g, chosenDisease, combo, otherOrder = 3, diseaseOrder = 0)

