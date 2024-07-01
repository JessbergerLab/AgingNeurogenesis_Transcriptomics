##== Setup the dataset ==##
######################################################################################################
# Load library
######################################################################################################
rm(list=ls())
library(dplyr)
library(MAST)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc)
library(biomaRt)
library(Seurat)
library(SingleCellExperiment)
library(ggrepel)
library(graphics)
library(Matrix)
library(lubridate)
library(glue)
library(gridExtra)
library(gdata)
library(slingshot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(gridExtra)
library(SummarizedExperiment)
library(DoubletFinder)
library(monocle)
######################################################################################################
# Load the dataset and add meatadata
######################################################################################################
# Read the data and create the Seurat Object
data.y <- Read10X_h5("./Data/Raw/Chromium/y.h5", use.names = TRUE, unique.features = TRUE)
data.m <- Read10X_h5("./Data/Raw/Chromium/m.h5", use.names = TRUE, unique.features = TRUE)
data.o1 <- Read10X_h5("./Data/Raw/Chromium/o1.h5", use.names = TRUE, unique.features = TRUE)
data.o2 <- Read10X_h5("./Data/Raw/Chromium/o2.h5", use.names = TRUE, unique.features = TRUE)
y <- CreateSeuratObject(counts = data.y, project = "y")
m <- CreateSeuratObject(counts = data.m, project = "m")
o1 <- CreateSeuratObject(counts = data.o1, project = "o1")
o2 <- CreateSeuratObject(counts = data.o2, project = "o2")
data.y <- NULL
data.m <- NULL
data.o1 <- NULL
data.o2 <- NULL
seurat_object <- merge(y, y = c(m, o1, o2), 
                       add.cell.ids = c("y", "m", "o1", "o2"), 
                       project = "AgingSpatialAtlas")
dim(seurat_object)
## [1] 32285 42478

# Add GEMgroup 
seurat_object@meta.data$GEMgroup[grepl("y_", rownames(seurat_object@meta.data))] <- "1"
seurat_object@meta.data$GEMgroup[grepl("m_", rownames(seurat_object@meta.data))] <- "2"
seurat_object@meta.data$GEMgroup[grepl("o1_", rownames(seurat_object@meta.data))] <- "3"
seurat_object@meta.data$GEMgroup[grepl("o2_", rownames(seurat_object@meta.data))] <- "4"
my_levels <- c("1", "2", "3", "4")
seurat_object@meta.data$GEMgroup <- factor(x = seurat_object@meta.data$GEMgroup, levels = my_levels)

# Add Age of samples 
seurat_object@meta.data$Age[grepl("y_", rownames(seurat_object@meta.data))] <- "Young"
seurat_object@meta.data$Age[grepl("m_", rownames(seurat_object@meta.data))] <- "Middle"
seurat_object@meta.data$Age[grepl("o1_|o2_", rownames(seurat_object@meta.data))] <- "Old"
my_levels <- c("Young", "Middle", "Old")
seurat_object@meta.data$Age <- factor(x = seurat_object@meta.data$Age, levels = my_levels)

# Add MT & Ribo content 
seurat_object[["percentMito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
seurat_object[["percentRibo"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rps|^Rpl")

# Add sex 
SexMatrix1 <- t(GetAssayData(seurat_object, assay = "RNA", slot = "counts"))
SexMatrix1 <- as.data.frame(SexMatrix1[ , "Xist"]) 
seurat_object[["sexScore"]] <- SexMatrix1
seurat_object[["Sex"]] <- ifelse(seurat_object[["sexScore"]] > 0, "Female", "Male") 




##== QC ==##
######################################################################################################
# QC
######################################################################################################
# Set thresholds for gene and mt content 
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "GEMgroup") 
VlnPlot(seurat_object, features = "percentMito", group.by = "GEMgroup", same.y.lims = TRUE, pt.size = 0.1) + ylim(0, 50)

# Save QC results as pictures
plot1 <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percentMito"), group.by = "GEMgroup", ncol = 3, pt.size = 0)
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "GEMgroup") 
plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percentMito", group.by = "GEMgroup") 

ggsave("QC1.png", plot = plot1, width = 15, height = 6)
ggsave("QC1.pdf", plot = plot1, width = 15, height = 6)
ggsave("QC2.png", plot = plot2)
ggsave("QC2.pdf", plot = plot2)
ggsave("QC3.png", plot = plot3)
ggsave("QC3.pdf", plot = plot3)

# Filter 
mitoHi <- 20 
nGeneLo <- 800
nCountLo <- 1000
nCountHi <- 30000
seurat_object <- subset(seurat_object, subset = nFeature_RNA > nGeneLo & nCount_RNA > nCountLo & nCount_RNA < nCountHi & percentMito < mitoHi)
dim(seurat_object)
## [1] 32285 39068





##== Doublet removal ==##
seurat_Y <- subset(seurat_object, idents = "Young")
#####################################################################################################
# Normalization and dimensionality reduction
#####################################################################################################
DefaultAssay(seurat_Y) <- "RNA"
seurat_Y <- NormalizeData(seurat_Y, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_Y <- RunUMAP(seurat_Y, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_Y, reduction = "umap")

#####################################################################################################
# DoubletFinder
#####################################################################################################
# pK identification
pc.num <- 1:30
sweep.res.list <- paramSweep_v3(seurat_Y, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic doublet proportion estimate
DoubletRate = 0.05  ## Assuming 5% doublet formation rate for 6000 cells
homotypic.prop <- modelHomotypic(seurat_Y$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(seurat_Y)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
seurat_Y <- doubletFinder_v3(seurat_Y, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                             nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# Visualization (results are stored in @meta.data)
DimPlot(seurat_Y, reduction = "umap", group.by = "DF.classifications_0.25_0.005_262")
seurat_Y <- SetIdent(seurat_Y, value = "DF.classifications_0.25_0.005_262")
seurat_Y %>% Idents() %>% table() # Check hw many cells in each cluster
## Singlet Doublet 
##    5516     262  

# Remove doublets from the data
seurat_Y <- subset(seurat_Y, idents = "Doublet", invert = TRUE)
seurat_Y <- SetIdent(seurat_Y, value = "seurat_clusters")

# Normalization and dimensionality reduction
seurat_Y <- NormalizeData(seurat_Y, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_Y <- RunUMAP(seurat_Y, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_Y, reduction = "umap")

# Save data
saveRDS(seurat_Y, "seurat_All_Y.rds")

seurat_M <- subset(seurat_object, idents = "Middle")
#####################################################################################################
# Normalization and dimensionality reduction
#####################################################################################################
DefaultAssay(seurat_M) <- "RNA"
seurat_M <- NormalizeData(seurat_M, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_M <- RunUMAP(seurat_M, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_M, reduction = "umap")

#####################################################################################################
# DoubletFinder
#####################################################################################################
# pK identification
pc.num <- 1:30
sweep.res.list <- paramSweep_v3(seurat_M, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic doublet proportion estimate
DoubletRate = 0.08  ## Assuming 8% doublet formation rate for 10000 cells
homotypic.prop <- modelHomotypic(seurat_M$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(seurat_M)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
seurat_M <- doubletFinder_v3(seurat_M, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                             nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# Visualization (results are stored in @meta.data)
DimPlot(seurat_M, reduction = "umap", group.by = "DF.classifications_0.25_0.005_753")
seurat_M <- SetIdent(seurat_M, value = "DF.classifications_0.25_0.005_753")
seurat_M %>% Idents() %>% table() # Check hw many cells in each cluster
## Singlet Doublet 
##    5359     419  

# Remove doublets from the data
seurat_M <- subset(seurat_M, idents = "Doublet", invert = TRUE)
seurat_M <- SetIdent(seurat_M, value = "seurat_clusters")

# Normalization and dimensionality reduction
seurat_M <- NormalizeData(seurat_M, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_M <- RunUMAP(seurat_M, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_M, reduction = "umap")

# Save data
saveRDS(seurat_M, "seurat_All_M.rds")

seurat_O <- subset(seurat_object, idents = "Old")
#####################################################################################################
# Normalization and dimensionality reduction
#####################################################################################################
DefaultAssay(seurat_O) <- "RNA"
seurat_O <- NormalizeData(seurat_O, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_O <- RunUMAP(seurat_O, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_O, reduction = "umap", label = T)

#####################################################################################################
# DoubletFinder
#####################################################################################################
# pK identification
pc.num <- 1:30
sweep.res.list <- paramSweep_v3(seurat_O, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Homotypic doublet proportion estimate
DoubletRate = 0.1  ## Assuming 10% doublet formation rate for over 10000 cells
homotypic.prop <- modelHomotypic(seurat_O$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(seurat_O)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
seurat_O <- doubletFinder_v3(seurat_O, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                             nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# Visualization (results are stored in @meta.data)
DimPlot(seurat_O, reduction = "umap", group.by = "DF.classifications_0.25_0.005_2053")
seurat_O <- SetIdent(seurat_O, value = "DF.classifications_0.25_0.005_2053")
seurat_O %>% Idents() %>% table() # Check hw many cells in each cluster
##  Singlet Doublet 
##    20892    2053  

# Remove doublets from the data
seurat_O <- subset(seurat_O, idents = "Doublet", invert = TRUE)
seurat_O <- SetIdent(seurat_O, value = "seurat_clusters")

DefaultAssay(seurat_O) <- "integrated"
seurat_O <- RunPCA(seurat_O, npcs = 50, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_O, reduction = "umap")

# Save data
saveRDS(seurat_O, "seurat_All_O.rds")

##== Data clustering ==##



sub_Young <- WhichCells(seurat_Y)
sub_Middle <- WhichCells(seurat_M)
sub_Old <- WhichCells(seurat_O)
seurat_object_counts <- subset(seurat_object, idents = c(sub_Young, sub_Middle,sub_Old))

# Assign Celltype (15 types)
seurat_object_counts$Doublet <- seurat_object_counts@active.ident
seurat_object_counts$Doublet <- "Yes"
Idents(seurat_object_counts) <- "Doublet"
Idents(seurat_object_counts, cells = sub_Young) <- "No" 
Idents(seurat_object_counts, cells = sub_Middle) <- "No" 
Idents(seurat_object_counts, cells = sub_Old) <- "No" 
seurat_object_counts <- subset(seurat_object_counts, idents = "No")

######################################################################################################
# Normalization & dimensionality reduction
######################################################################################################
DefaultAssay(seurat_object_counts) <- "RNA"
seurat_object.integrated <- NormalizeData(seurat_object_counts, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
DimPlot(seurat_object.integrated, reduction = "umap")

######################################################################################################
# Assign cell cycle genes
######################################################################################################
CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(seurat_object.integrated))
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(seurat_object.integrated))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(seurat_object.integrated))
seurat_object.integrated <- CellCycleScoring(object = seurat_object.integrated, g2m.features = g2m_genes, s.features = s_genes)
FeaturePlot(seurat_object.integrated, features = c("S.Score", "G2M.Score"), 
            pt.size = 0.5, cols = c(rgb(0.7, 0.7, 0.7, alpha = 0.2), "darkred"), min.cutoff = "q10", max.cutoff = "q90")




