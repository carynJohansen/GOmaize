# Introduction and About ------------

# Libraries and Set Up --------------

# bioconductor packages
# topGO
# GO.db

# reduce the maize.db to only the maize data

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

# Main ------------------------------

# make an analysis folder

# AgriGO provides it's maize genes in v3 gene name format. 
# Not sure what to do about this right now, except convert my v4 genes back to v3 names

# ddplyr filter
# collect names of genes lost in this transition

# merge/union of v3 gene list and db
# collect genes lost in this step

# total percent of genes from the original list now with GO terms
# is this hig enough to justify this approach?

# create list of genes for GO term analysis

# create the topGOdata object

# classical enrichment analysis - runTest()
#save this result

# Kolmogorov-Smirnov test, classic method
# save this result

# Kolmogorov-Smirnov test, elim method
# save this result

# use GenTable to get a dataframe containing the top GO nodes

# Output ----------------------------

# I want to save:
# - genes that got lost along the way. How many and what are they? During which step were they lost?
# - GO for genes remaining, the resutls of the GenTable perhaps