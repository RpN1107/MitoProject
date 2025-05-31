
# --------------------------------MITOPROJECT-----------------------------------
# Rithwik Nambiar-----------------------------------------------OMICS Assignment 
# 20244013-----------------------------------------------------------------BIODS

# ------------------------------------------------------------------------------
# Objective: 
# To perform GO and reactome analysis on the top 10 features extracted after ML

# Data Used: top_mito_genes_for_GO.csv

# ------------------------------------------------------------------------------

# Step 0 : Setup working directory
setwd("D:/BioDS/Sem2/OMICS_Analysis/Project/3_Stimulated and frozen PBMCs")
set.seed(42)

# ------------------------------------------------------------------------------

# Step 1: Installing and loading the necessary packages

#install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# ------------------------------------------------------------------------------
# Step 2: Preparing the gene list

# Ensure gene names are symbols (e.g., "COX6C", "NDUFA9", etc.)
df <- read.csv("top_mito_genes_for_GO.csv")
top10_genes <- df[,1]

# Convert gene symbols to Entrez IDs (required for most GO/Reactome tools)
entrez_ids <- bitr(top10_genes, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ------------------------------------------------------------------------------
# Step 3: GO for Biological Processes

ego_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

dotplot(ego_bp) + ggtitle("Top 10 features - Biological Process")

# ------------------------------------------------------------------------------
# Step 4: GO for Molecular Function
ego_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

dotplot(ego_mf) + ggtitle("Top 10 features - Molecular Function")

# ------------------------------------------------------------------------------
# Step 5: GO for Cellular Component
ego_cc <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

dotplot(ego_cc) + ggtitle("Top 10 features - Cellular Component")

# ------------------------------------------------------------------------------
# Step 6: Enrichment for Reactome Pathway
reactome_results <- enrichPathway(gene = entrez_ids$ENTREZID,
                                  organism = "human",
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)

dotplot(reactome_results) + ggtitle("Top 10 features - Reactome Pathway")

# ------------------------------------------------------------------------------
