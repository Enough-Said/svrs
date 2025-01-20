#
# Clean the drug combinations data and convert the IDs
#

library(dplyr)
library(tidyr)
library(jsonlite)

knownCombs <- read.csv("rawdata/knownPairwiseDrugCombs.csv", header = TRUE)
mapID <- read.csv("rawdata/drug-mappings.tsv", header = TRUE, sep = "\t")
mapID <- mapID[c("drugbankId", "chembl_id")] %>% filter(chembl_id != "null")

# Convert drugbank IDs to CHEMBL IDs for consistency with the network
knownCombs$interactionA <- lapply(knownCombs[, 1], function(x) mapID[mapID$drugbankId == x, "chembl_id"])
knownCombs$interactionB <- lapply(knownCombs[, 2], function(x) mapID[mapID$drugbankId == x, "chembl_id"])

# Expand one-to-many mappings to include all possible useful maps
knownCombs <- knownCombs %>% 
    unnest(cols = c(interactionA, interactionB)) %>% 
    as.data.frame()
knownCombs <- knownCombs[c("interactionA", "interactionB")]

write(toJSON(knownCombs), "clean/verifiedCombinations.json")
