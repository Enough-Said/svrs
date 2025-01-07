#
# A pipeline displaying all the scripts and functions together to find
# drug combinations
#

source("analysisScripts/filterNodes.r")
source("analysisScripts/filterEdges.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")

# tissue <- fromJSON("clean/tissue.json")
# cellLine <- read.delim("clean/cellLine.tsv")
# cellType <- fromJSON("clean/cellType.json")
# subcell <- fromJSON("clean/subcell.json")

################################################################################

chosenDisease <- "ALZHEIMER DISEASE 2"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.tissue(g)

drugDist <- getDistFromDisease(g, chosenDisease, 1)
drugPairs <- findDrugPairs(g, chosenDisease, 1, drugDist)
drugCombs <- findDrugCombinations(g, 
    chosenDisease, 1, precompPairs = drugPairs, 
    maxSize = 2)

print(drugCombs[[1]][[2]])

################################################################################

chosenDisease <- "Alzheimer's Disease"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.type(g)
newg <- v.filter.tissue(g, "hippocampal formation", minNTPM = 1)

drugPairs <- findDrugPairs(newg, chosenDisease, 0)
drugCombs <- findDrugCombinations(newg, chosenDisease, 0, precompPairs = drugPairs)

print(drugCombs)
