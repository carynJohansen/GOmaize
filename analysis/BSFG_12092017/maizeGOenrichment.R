# Introduction and About ------------

# Libraries and Set Up --------------

setwd("/Users/caryn/Box Sync/Projects/GOmaize/")

# bioconductor packages
# topGO
# GO.db

# reduce the maize.db to only the maize data

library(dplyr)
library(topGO)

# Data ------------------------------

# loaded data will include:
# - list of genes of interest, this is gene IDs and perhaps a p-value/cor value/something based on previous analysis
# - maize database from AgriGO
# - converstion table for v3 to v4

# For topGO to work, I need to create a topGOdata object. this object needs:
# - list of gene identifiers and gene-wise scores (optional) such as the p-values from a DE analysis
# - a mapping btwn gene IDs and GO terms. this is the annotation db downloaded from AgriGO
# - GO heirarichal strucutre, from the GO.db package
# - a list of genes to define the gene universe
# - an annotation function, here I will be using annFUN.gene2GO
# - gene2GO annotation function has required arguments

genesv3 <- readRDS("analysis/BSFG_12092017/modules_09122017.rds")

Ngenes <- n_distinct(genesv3[,1])

# these genes are already v3 or v2
#v4tov3 <- read.table("db/MaizeGDB_v3_v4.genes.csv", header=T, sep=",")

maize.db <- read.table("db/maize.agrigo.f.db", sep="\t", header=F)
#maize.db <- dplyr::filter(maize.db, source=="maize")

agrigo <- read.table("db/maize_v3.agrigo.f.db", header=F, sep="\t")
argot <- read.table("db/maize_v3.argot2.f.db", header=F, sep="\t")

# Main ------------------------------

# make an analysis folder

analysisDir <- "BSFG_12092017"
#dir.create(file.path(paste("analysis",analysisDir, sep="/")))

# AgriGO provides it's maize genes in v3 gene name format. 
# Not sure what to do about this right now, except convert my v4 genes back to v3 names

# ddplyr filter
# collect names of genes lost in this transition

#genesv4tov3 <- v4tov3[v4tov3$v4.gene.ID..if.present. %in% genesv4$genes,]

#genesv3 <- dplyr::select(genesv4tov3, v3.gene.ID,  v4.gene.ID..if.present.)
#genesv3 <- merge(genesv3, genesv4, by = "v4.gene.ID..if.present.")#, by.y = "gene")
#genesv3$v4.gene.ID..if.present. <- NULL #remove the v4 id column

#lostGenes_v4tov3 <- genesv4[!c(genesv4$v4.gene.ID..if.present. %in% v4tov3_sm$v4.gene.ID..if.present.),]

# merge/union of v3 gene list and db
# collect genes lost in this step

genes.in.db <- maize.db[maize.db$V1 %in% genesv3$genes,]
n.genes.in.db <- n_distinct(genes.in.db$V1)

# total percent of genes from the original list now with GO terms
# is this hig enough to justify this approach?
percent_N_go_analysis <- n.genes.in.db/Ngenes

# create list of genes for GO term analysis

#maize.db$source <- NULL
#write.table(maize.db, file="db/maize.agrigo.f.db", sep="\t", quote = FALSE, row.names = F, col.names = F)
#maize.db <- read.table("db/maize.agrigo.f.db", header= F, sep="\t")

maize.gene2GO <- readMappings(file = "db/maize.agrigo.all.f.db")

agrigo.gene2GO <- readMappings(file= "db/maize_v3.agrigo.f.db")

argot.gene2GO <- readMappings(file="db/maize_v3.argot2.f.db")

#maize.gene2GO <- maize.gene2GO[names(maize.gene2GO) %in% genesv3$v3.gene.ID]
#genesv3 <- genesv3[genesv3$v3.gene.ID %in% names(maize.gene2GO),]

# Establish the genes of interest

# Function to select genes for GO enrichment, based on the threshold value, in the `value` column
goi <- function(thresh) {
  return((thresh < -0.2 ) |  (thresh > 0.2)) 
}

#For this analysis, I need to do a GO term analysis for each module

n_modules <- n_distinct(genesv3$factor)

module_GO_1 <- list()

for (m in 1:n_modules) {
  mlist <- list()
  genes <- filter(genesv3, factor == m)
  
  #full join to get the factor values associated with the right genes in the db
  maize.db.m <- full_join(maize.db, genes, by = c("V1" = "genes"))
  
  # clean db, get rid of genes not in db
  lostGenes <- maize.db.m[is.na(maize.db.m$V2),]
  maize.db.m <- maize.db.m[!is.na(maize.db.m$V2),]
  
  #replace the NA's in the value column with zeros
  maize.db.m[is.na(maize.db.m)] <- 0
  
  # get geneList for topGO
  geneList <- maize.db.m$value
  names(geneList) <- maize.db.m$V1
  
  #build the GO data term
  GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, 
                geneSel = goi,
                nodeSize = 5,
                annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_MF, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_MF, algorithm = "elim", statistic = "ks")
  allRes_MF <- GenTable(GOdata_MF, classicFisher = result_fisher, 
                     classicKS = result_KS_classic,
                     elimKS = result_KS_elim,
                     orderBy = "classicFisher",
                     ranksOf = "classicFisher",
                     topNodes = 100)
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                   geneSel = goi,
                   nodeSize = 5,
                   annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_BP, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_BP, algorithm = "elim", statistic = "ks")
  allRes_BP <- GenTable(GOdata_BP, classicFisher = result_fisher, 
                        classicKS = result_KS_classic,
                        elimKS = result_KS_elim,
                        orderBy = "classicFisher",
                        ranksOf = "classicFisher",
                        topNodes = 100)
  
  # save results
  mlist <- list(genes, lostGenes, maize.db.m, geneList, GOdata_MF, allRes_MF, GOdata_BP, allRes_BP)
  module_GO_1[[m]] <- mlist
}

# for the agrigo from GAMER

module_GO_2 <- list()

for (m in 1:n_modules) {
  mlist <- list()
  genes <- filter(genesv3, factor == m)
  
  #full join to get the factor values associated with the right genes in the db
  agrigo.m <- full_join(agrigo, genes, by = c("V1" = "genes"))
  
  # clean db, get rid of genes not in db
  lostGenes <- agrigo.m[is.na(agrigo.m$V2),]
  agrigo.m <- agrigo.m[!is.na(agrigo.m$V2),]
  
  #replace the NA's in the value column with zeros
  agrigo.m[is.na(agrigo.m)] <- 0
  
  # get geneList for topGO
  geneList <- agrigo.m$value
  names(geneList) <- agrigo.m$V1
  
  #build the GO data term
  GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, 
                   geneSel = goi,
                   nodeSize = 5,
                   annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_MF, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_MF, algorithm = "elim", statistic = "ks")
  allRes_MF <- GenTable(GOdata_MF, classicFisher = result_fisher, 
                        classicKS = result_KS_classic,
                        elimKS = result_KS_elim,
                        orderBy = "classicFisher",
                        ranksOf = "classicFisher",
                        topNodes = 100)
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                   geneSel = goi,
                   nodeSize = 5,
                   annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_BP, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_BP, algorithm = "elim", statistic = "ks")
  allRes_BP <- GenTable(GOdata_BP, classicFisher = result_fisher, 
                        classicKS = result_KS_classic,
                        elimKS = result_KS_elim,
                        orderBy = "classicFisher",
                        ranksOf = "classicFisher",
                        topNodes = 100)
  
  # save results
  mlist <- list(genes, lostGenes, agrigo.m, geneList, GOdata_MF, allRes_MF, GOdata_BP, allRes_BP)
  module_GO_2[[m]] <- mlist
}

# Using the GAMER annotations

module_GO_1 <- list()

for (m in 1:n_modules) {
  mlist <- list()
  genes <- filter(genesv3, factor == m)
  
  #full join to get the factor values associated with the right genes in the db
  maize.db.m <- full_join(maize.db, genes, by = c("V1" = "genes"))
  
  # clean db, get rid of genes not in db
  lostGenes <- maize.db.m[is.na(maize.db.m$V2),]
  maize.db.m <- maize.db.m[!is.na(maize.db.m$V2),]
  
  #replace the NA's in the value column with zeros
  maize.db.m[is.na(maize.db.m)] <- 0
  
  # get geneList for topGO
  geneList <- maize.db.m$value
  names(geneList) <- maize.db.m$V1
  
  #build the GO data term
  GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, 
                   geneSel = goi,
                   nodeSize = 5,
                   annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_MF, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_MF, algorithm = "elim", statistic = "ks")
  allRes_MF <- GenTable(GOdata_MF, classicFisher = result_fisher, 
                        classicKS = result_KS_classic,
                        elimKS = result_KS_elim,
                        orderBy = "classicFisher",
                        ranksOf = "classicFisher",
                        topNodes = 100)
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                   geneSel = goi,
                   nodeSize = 5,
                   annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)
  
  # tests
  result_fisher <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
  result_KS_classic <- runTest(GOdata_BP, algorithm = "classic", statistic = "ks")
  result_KS_elim <- runTest(GOdata_BP, algorithm = "elim", statistic = "ks")
  allRes_BP <- GenTable(GOdata_BP, classicFisher = result_fisher, 
                        classicKS = result_KS_classic,
                        elimKS = result_KS_elim,
                        orderBy = "classicFisher",
                        ranksOf = "classicFisher",
                        topNodes = 100)
  
  # save results
  mlist <- list(genes, lostGenes, maize.db.m, geneList, GOdata_MF, allRes_MF, GOdata_BP, allRes_BP)
  module_GO_1[[m]] <- mlist
}

# Output ----------------------------

# I want to save:
# - genes that got lost along the way. How many and what are they? During which step were they lost?
# - GO for genes remaining, the resutls of the GenTable perhaps

save(module_GO, file = paste(paste("analysis/",analysisDir,sep=""), paste(analysisDir, ".RData", sep=""), sep="/"))
