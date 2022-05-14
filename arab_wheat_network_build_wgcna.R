library(WGCNA)
​
#arabidopsis
arab_net = blockwiseModules(data,
                            networkType = "signed",
                            power = 18,
                            TOMType = "signed",
                            minModuleSize = 30,
                            corType = "bicor",
                            maxPOutliers = 0.05,
                            reassignThreshold = 0,
                            mergeCutHeight = 0.15,
                            numericLabels = TRUE,
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = paste(out_filename, ".Rdata", sep=""),
                            verbose = 3)
​
#################################################################################################################
#################################################################################################################
#wheat 
wheat_net = blockwiseModules(data, 
                       networkType = "signed",
                       power = 18,
                       TOMType = "signed",
                       minModuleSize = 30,
                       corType = "bicor",
                       maxPOutliers = 0.05,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.15,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste(out_filename, ".Rdata", sep=""),
                       verbose = 3)
#################################################################################################################
#Call an automatic merging function
merge = mergeCloseModules(data, module_colours_for_genes, MEs = MEs, cutHeight = 0.25, verbose = 5, iterate = F)
​
#The merged module colors
newmergedColors = merge$colors
​
# Eigengenes of the new merged modules:
newMEs = merge$newMEs
​
#Rename to moduleColors
mrg_moduleColors = newmergedColors
##################################################################################################################
