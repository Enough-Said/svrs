#
# A pipeline displaying all the scripts and functions together to find
# drug combinations
#

source("scripts/filterNodes.r")
source("scripts/filterEdges.r")
source("scripts/findDistance.r")
source("scripts/findDrugs.r")

tissue <- fromJSON("clean/tissue.json")
cellLine <- read.delim("clean/cellLine.tsv")
cellType <- fromJSON("clean/cellType.json")
subcell <- fromJSON("clean/subcell.json")

################################################################################

chosenDisease <- "ALZHEIMER DISEASE 2"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.tissue(g)

drugDist <- getDistFromDisease(g, chosenDisease, 0)
drugPairs <- findDrugPairs(g, chosenDisease, 0, drugDist)
drugCombs <- findDrugCombinations(g, chosenDisease, 0, precompPairs = drugPairs)

print(drugCombs)

################################################################################

chosenDisease <- "Alzheimer's Disease"

g <- readRDS("clean/baseGraph.rds")
g <- e.filter.subcell(g)
g <- e.filter.type(g)
newg <- v.filter.tissue(g, "hippocampal formation", minNTPM = 1)

drugPairs <- findDrugPairs(newg, chosenDisease, 0)
drugCombs <- findDrugCombinations(newg, chosenDisease, 0, precompPairs = drugPairs)

print(drugCombs)
