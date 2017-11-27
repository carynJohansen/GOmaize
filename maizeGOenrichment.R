# Introduction and About ------------

# Libraries and Set Up --------------

setwd("/Users/caryn/Box Sync/Projects/GOmaize/")

# bioconductor packages
# topGO
# GO.db

# reduce the maize.db to only the maize data

library(dplyr)

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

genesv4 <- readRDS("analysis/BSFG_12092017/modules_09122017.rds")

Ngenes <- n_distinct(genesv4[,1])

v4tov3 <- read.table("db/MaizeGDB_v3_v4.genes.csv", header=T, sep=",")

maize.db <- read.table("db/maize.agrigo.f.db", sep="\t", header=T)
#maize.db <- dplyr::filter(maize.db, source=="maize")

# Main ------------------------------

# make an analysis folder

analysisDir <- "B73_TIL01_DE"
dir.create(file.path(paste("analysis",analysisDir, sep="/")))

# AgriGO provides it's maize genes in v3 gene name format. 
# Not sure what to do about this right now, except convert my v4 genes back to v3 names

# ddplyr filter
# collect names of genes lost in this transition

colnames(genesv4) <- c("v4.gene.ID..if.present.", "Adj.P.Value")
genesv4tov3 <- v4tov3[v4tov3$v4.gene.ID..if.present. %in% genesv4[,1],]

genesv3 <- dplyr::select(genesv4tov3, v3.gene.ID,  v4.gene.ID..if.present.)
genesv3 <- merge(genesv3, genesv4, by = "v4.gene.ID..if.present.")#, by.y = "gene")
genesv3$v4.gene.ID..if.present. <- NULL #remove the v4 id column

#lostGenes_v4tov3 <- genesv4[!c(genesv4$v4.gene.ID..if.present. %in% v4tov3_sm$v4.gene.ID..if.present.),]

# merge/union of v3 gene list and db
# collect genes lost in this step

genes.in.db <- maize.db[maize.db$V1 %in% genesv3$v3.gene.ID,]
n.genes.in.db <- n_distinct(genes.in.db$V1)

# total percent of genes from the original list now with GO terms
# is this hig enough to justify this approach?
percent_N_go_analysis <- n.genes.in.db/Ngenes

# create list of genes for GO term analysis

#maize.db$source <- NULL
#write.table(maize.db, file="db/maize.agrigo.f.db", sep="\t", quote = FALSE, row.names = F, col.names = F)
maize.db <- read.table("db/maize.agrigo.f.db", header= F, sep="\t")
maize.gene2GO <- readMappings(file = "db/maize.agrigo.f.db")

#maize.gene2GO <- maize.gene2GO[names(maize.gene2GO) %in% genesv3$v3.gene.ID]
#genesv3 <- genesv3[genesv3$v3.gene.ID %in% names(maize.gene2GO),]

# Establish the genes of interest
genesv3_list <- maize.db$V1 %in% genesv3$v3.gene.ID
maize.db$pvalue <- 1

for (i in 1:nrow(genesv3)) {
  for (j in 1:nrow(maize.db)) {
    if (as.character(genesv3$v3.gene.ID[i]) == as.character(maize.db$V1[j])) {
      maize.db$pvalue[j] <- genesv3$Adj.P.Value[i]
    }
  }
}


geneList <- maize.db$pvalue
names(geneList) <- maize.db$V1
#geneList <- genesv3$adj.P.Val.x
#names(geneList) <- genesv3$v3.gene.ID

# Select Gene function

DiffGenes <- function(pvalue) {
  return(pvalue < 0.01)
}

goi <- function(pvalue) {
  return(pvalue < 0.01)
}

x <- DiffGenes(geneList)
sum(x)
# create the topGOdata object

x <- goi(geneList)
sum(x)

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, 
              geneSel = goi,
              nodeSize = 5,
              annotationFun = annFUN.gene2GO, gene2GO = maize.gene2GO)

# classical enrichment analysis - runTest()
#save this result

result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Kolmogorov-Smirnov test, classic method
# save this result

result_KS_classic <- runTest(GOdata, algorithm = "classic", statistic = "ks")

# Kolmogorov-Smirnov test, elim method
# save this result

result_KS_elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# use GenTable to get a dataframe containing the top GO nodes

(allRes <- GenTable(GOdata, classicFisher = result_fisher, 
                   classicKS = result_KS_classic,
                   elimKS = result_KS_elim,
                   orderBy = "elimKS",
                   ranksOf = "classicFisher",
                   topNodes = 30))

# Output ----------------------------

# I want to save:
# - genes that got lost along the way. How many and what are they? During which step were they lost?
# - GO for genes remaining, the resutls of the GenTable perhaps

save(maize.db, genesv3, GOdata, allRes, file = paste(analysisDir, paste(analysisDir, ".RData", sep=""), sep="/"))
