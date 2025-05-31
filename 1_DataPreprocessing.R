
# --------------------------------MITOPROJECT-----------------------------------
# Rithwik Nambiar-----------------------------------------------OMICS Assignment 
# 20244013-----------------------------------------------------------------BIODS

# ------------------------------------------------------------------------------
# Objective: 
# 1. To process the scRNA-seq data for ML 
# 2. Compare the expression patterns of nuclear and mitochondrial encoded genes
#    in control vs LPS-treated conditions

# Data Used: PMID: 37414801 Series GSE226488	
# Samples: GSM7077866 (lps) GSM7077867 (control)

# ------------------------------------------------------------------------------

# Step 0 : Setup working directory
setwd("D:/BioDS/Sem2/OMICS_Analysis/Project/MitoProject_Rithwik/")

# ------------------------------------------------------------------------------

# Step 1: Installing and loading the necessary packages

#install.packages("Seurat")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("readxl")
#install.packages("BiocManager")
#BiocManager::install("celldex")
#BiocManager::install("SingleR")

library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(celldex)
library(SingleR)
# ------------------------------------------------------------------------------

# Step 2: Preparing the dataset

# 2.1 Reading the data
control_data <- Read10X(data.dir="./Input/control/")
lps_data <- Read10X(data.dir="./Input/lps/")

# 2.2 Creating the Seurat Objects
control_seurat <- CreateSeuratObject(counts = control_data, project = "control")
lps_seurat <- CreateSeuratObject(counts = lps_data, project = "LPS")

# 2.3 Adding the condition labels
control_seurat$condition <- "control"
lps_seurat$condition <- "LPS"

# 2.4 Merging Both datasets
combined_seurat <- merge(control_seurat, y = lps_seurat, add.cell.ids = c("control", "LPS"))

# ------------------------------------------------------------------------------

# Step 3: Quality Control and Filtering 

# 3.1 Add percent.mt
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")

# 3.2 Basic QC 
VlnPlot(combined_seurat, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 3.3 Filtering
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 1000 & 
                  nFeature_RNA < 7500 & percent.mt < 30 )

# ------------------------------------------------------------------------------

# Step 4: Normalisation and Scaling

# 4.1 Normalize the data using LogNormalize
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# 4.2 Identify the top 2,000 variable features
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)

# 4.3 Scale the data
combined_seurat <- ScaleData(combined_seurat, features = rownames(combined_seurat))

# ------------------------------------------------------------------------------
# Step 5: Dimensionality Reduction

# 5.1 Dimensionality Reduction
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))

# ------------------------------------------------------------------------------
# Step 6: Filtering the Dataset to only have the mitochondrial genes

# 6.1 Set correct default assay
DefaultAssay(combined_seurat) <- "RNA"

# 6.2 Access log-normalized expression correctly from 'data' layer
control_data <- combined_seurat[["RNA"]]@layers[["data.control"]]
lps_data   <- combined_seurat[["RNA"]]@layers[["data.LPS"]]

# 6.3 Assigning Barcodes properly
# control
control_barcodes <- readLines("control/barcodes.tsv.gz")
control_features <- read.delim("control/features.tsv.gz", header = FALSE)

# LPS
lps_barcodes <- readLines("lps/barcodes.tsv.gz")
lps_features <- read.delim("lps/features.tsv.gz", header = FALSE)

control_barcodes <- control_barcodes[1:ncol(control_data)]
lps_barcodes   <- lps_barcodes[1:ncol(lps_data)]

# Assign to control_data
rownames(control_data) <- control_features$V2  # or V1 if V2 is missing
colnames(control_data) <- paste0("control_", control_barcodes)

# Assign to lps_data
rownames(lps_data) <- lps_features$V2
colnames(lps_data) <- paste0("LPS_", lps_barcodes)

# 6.4 Get mitochondrial and nuclear encoded mitochondrial genes

#Mitochondrial genes
mt_genes <- grep("^MT-", rownames(combined_seurat), value = TRUE)

#Genes taken from MitoCart3.0
mitocarta <- read_excel("Human.MitoCarta3.0.xls", sheet = 2)
colnames(mitocarta)
mitocarta_genes <- unique(mitocarta$Symbol) 

all_mito_genes <- unique(c(mt_genes, mitocarta_genes)) #Combining both gene_sets

#Getting OXPHOS genes
mitocarta_OXPHOS <- read_excel("Human.MitoCarta3.0.xls", sheet = 4) 

# Check column names
colnames(mitocarta_OXPHOS)

#Extracting OXPHOS genes
oxphos_genes <- unique(mitocarta_OXPHOS$Genes[grepl("OXPHOS", mitocarta_OXPHOS$MitoPathway)])
oxphos_genes <- trimws(unlist(strsplit(paste(oxphos_genes, collapse = ","), ",")))

# Filtering genes that exist in the seurat object
oxphos_genes_in_seurat <- intersect(rownames(combined_seurat), oxphos_genes)
# ------------------------------------------------------------------------------
# step 7: Exporting the mito data

# 7.1 Filtering the genes that exist in the seurat object
mito_genes_in_seurat <- intersect(rownames(combined_seurat), all_mito_genes)

# 7.2 Subsetting the seurat object to only include mito genes
seurat_mito_only <- subset(combined_seurat, features = mito_genes_in_seurat)

# 7.3 Visualising the data
VlnPlot(seurat_mito_only, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 7.4 Extracting the expression data
mito_data <- FetchData(seurat_mito_only, vars = mito_genes_in_seurat)
#Adding condition column (to be used as label)
mito_data$condition <- seurat_mito_only@meta.data[rownames(mito_data), "condition"]

# 7.5 Writing the excel file
write.csv(mito_data, "Mito_expression_data.csv")

# ------------------------------------------------------------------------------
# step 8: Clustering the Data

# 8.1 Run PCA and Find Neighbours and clusters
seurat_mito_only <- RunPCA(seurat_mito_only, features = VariableFeatures(seurat_mito_only))
seurat_mito_only <- FindNeighbors(seurat_mito_only, dims = 1:10)
seurat_mito_only <- FindClusters(seurat_mito_only, resolution = 0.2) 

# 8.2 Visulaising the plot
DimPlot(seurat_mito_only, reduction = "pca", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("PCA Clustering")

# 8.3 Only comparing clusters 1 an 3
seurat_mito_only <- JoinLayers(seurat_mito_only)
seurat_mito_subset <- subset(seurat_mito_only, idents = c(1, 3))

# 8.4 Finding the top 10 markers in clusters 1 and 3
markers <- FindAllMarkers(seurat_mito_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# 8.5 Visualising
DoHeatmap(seurat_mito_subset, features = top10$gene)

# ------------------------------------------------------------------------------
# Step 9: Comparing nuclear and mitochondrial encoded genes

# 9.1 Get the genes in clusters 1 and 3
genes_in_data <- rownames(seurat_mito_subset)

# 9.2 Select the oxphos and mt genes in those clusters
oxphos_filtered <- intersect(oxphos_genes, genes_in_data)
mt_filtered <- intersect(mt_genes, genes_in_data)

# 9.3 Extract the most variable genes between conditions
Idents(seurat_mito_subset) <- "condition"

#OXPHOS
oxphos_de <- FindMarkers(seurat_mito_subset, ident.1 = "control", ident.2 = "LPS",
                         features = oxphos_filtered, only.pos = FALSE,
                         min.pct = 0.1, logfc.threshold = 0)

mt_de <- FindMarkers(seurat_mito_subset, ident.1 = "control", ident.2 = "LPS",
                     features = mt_filtered, only.pos = FALSE,
                     min.pct = 0.1, logfc.threshold = 0)

# 9.4 Select the most variable OXPHOS and MT genes
#OXPHOS
top_oxphos <- oxphos_de %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5)

#MT
top_mt <- mt_de %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5)

# 9.5 Visualising these genes
top_genes <- unique(c(rownames(top_oxphos), rownames(top_mt)))

VlnPlot(seurat_mito_subset, features = top_genes, 
        group.by = "condition", 
        pt.size = 0.1)

# ------------------------------------------------------------------------------
# Step 10: Identify which PBMCs are present in clusters 1 and 3 

# Convert Seurat object to SingleCellExperiment
expr_matrix <- GetAssayData(seurat_mito_subset, layer = "data")

# Use Human Primary Cell Atlas as reference
ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = expr_matrix, ref = ref, labels = ref$label.main)

# Add labels back to Seurat object
seurat_mito_subset$SingleR_labels <- pred$labels

# Visualize annotations
DimPlot(seurat_mito_subset, group.by = "SingleR_labels", label = TRUE, reduction = "pca") +
  ggtitle("Cell Type Annotation with SingleR")

# ------------------------------------------------------------------------------