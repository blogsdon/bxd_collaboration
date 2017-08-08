#### pull rosmap modules
synapseClient::synapseLogin()
moduleDefs <- synapseClient::synTableQuery("select * from syn10226075")@values

#### load data from dropbox
mouseModules <- load('~/Dropbox/AD_BXD/AllSamples/BXDOnly_14months_allsamples.RData')

#### reorganize mouse modules
mouseModuleDefn <- data.frame(GeneID=colnames(datmonths14),
                              Module=moduleColors,
                              stringsAsFactors=FALSE)


#### get mouse orthologs
moduleDefs2 <- utilityFunctions::getMouseOrthologsFromHumanIds(moduleDefs$GeneID)

#### reduce down to common sets
combinedManifest <- dplyr::select(moduleDefs,GeneID,ModuleNameFull)
combinedManifest <- dplyr::left_join(moduleDefs2,combinedManifest,by=c('ensembl_gene_id'='GeneID'))
combinedManifest <- dplyr::left_join(combinedManifest,mouseModuleDefn,by=c('mmusculus_homolog_ensembl_gene'='GeneID'))

#### remove dups and missing
which(is.na(combinedManifest$Module))
combinedManifest <- combinedManifest[which(!is.na(combinedManifest$Module)),]

#### listify
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
rosmapModulesEnsg <- lapply(unique(combinedManifest$ModuleNameFull),listify,combinedManifest$ensembl_gene_id,combinedManifest$ModuleNameFull)
names(rosmapModulesEnsg) <- unique(combinedManifest$ModuleNameFull)
rosmapModulesEnsg2 <- lapply(rosmapModulesEnsg,unique)


mouseModulesEnsg <- lapply(unique(combinedManifest$Module),listify,combinedManifest$ensembl_gene_id,combinedManifest$Module)
names(mouseModulesEnsg) <- unique(combinedManifest$Module)
mouseModulesEnsg2 <- lapply(mouseModulesEnsg,unique)

human_mouse_or <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,
                                      rosmapModulesEnsg2,
                                      mouseModulesEnsg2,
                                      unique(unlist(rosmapModulesEnsg2)))

human_mouse_pval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                      rosmapModulesEnsg2,
                                      mouseModulesEnsg2,
                                      unique(unlist(rosmapModulesEnsg2)))




rosmapModules <- lapply(unique(combinedManifest$ModuleNameFull),listify,combinedManifest$external_gene_name,combinedManifest$ModuleNameFull)
names(rosmapModules) <- unique(combinedManifest$ModuleNameFull)
rosmapModules2 <- lapply(rosmapModules,unique)


mouseModules <- lapply(unique(combinedManifest$Module),listify,combinedManifest$external_gene_name,combinedManifest$Module)
names(mouseModules) <- unique(combinedManifest$Module)
mouseModules2 <- lapply(mouseModules,unique)
mouseModules2 <- lapply(mouseModules2,toupper)


genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]

cellList <- GeneSets$Cell_Markers
cellList$ADGWAS <- adList

human_ad_or <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,
                                                rosmapModules2,
                                                cellList,
                                                unique(unlist(rosmapModules2)))

human_ad_pval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                                  rosmapModules2,
                                                  cellList,
                                                  unique(unlist(rosmapModules2)))

mouse_ad_or <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,
                                             mouseModules2,
                                             cellList,
                                             unique(unlist(rosmapModules2)))

mouse_ad_pval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                               mouseModules2,
                                               cellList,
                                               unique(unlist(rosmapModules2)))

