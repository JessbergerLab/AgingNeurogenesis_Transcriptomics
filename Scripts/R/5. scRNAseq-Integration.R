######################################################################################################
# Load the data
######################################################################################################
seurat1 <- readRDS("./Data/SeuratObject/Chromium/Total/Annotation/Final_Annotation.rds")
seurat2 <- readRDS("./Data/SeuratObject/Chromium/Integration/Hochgerner_Dataset_C.rds")

######################################################################################################
# Prepare the metadata
######################################################################################################
# Adjust GEMgroup infomation
seurat1$BBB <- seurat1$GEMgroup
seurat1$BBB <- gsub("1", "a1", seurat1$BBB)
seurat1$BBB <- gsub("2", "a2", seurat1$BBB)
seurat1$BBB <- gsub("3", "a3", seurat1$BBB)
seurat1$BBB <- gsub("4", "a4", seurat1$BBB)
seurat1$GEMgroup <- seurat1$BBB
seurat1$BBB <- NULL

Idents(seurat2) <- seurat2$GEMgroup
seurat2$BBB <- seurat2$GEMgroup
seurat2$BBB <- gsub("1", "b1", seurat2$BBB)
seurat2$BBB <- gsub("2", "b2", seurat2$BBB)
seurat2$BBB <- gsub("3", "b3", seurat2$BBB)
seurat2$BBB <- gsub("4", "b4", seurat2$BBB)
seurat2$BBB <- gsub("5", "b5", seurat2$BBB)
seurat2$BBB <- gsub("6", "b6", seurat2$BBB)
seurat2$BBB <- gsub("7", "b7", seurat2$BBB)
seurat2$BBB <- gsub("8", "b8", seurat2$BBB)
seurat2$BBB <- gsub("9", "b9", seurat2$BBB)
seurat2$BBB <- gsub("10", "b10", seurat2$BBB)
seurat2$BBB <- gsub("11", "b11", seurat2$BBB)
seurat2$GEMgroup <- seurat2$BBB
seurat2$BBB <- NULL

# Extract metadata information from our DG
metadata1 <- data.frame(seurat1@meta.data)
Cellid <- rownames(metadata1)
metadata1 <- cbind(metadata1, new_col = Cellid)
colnames(metadata1)[19] <- "Cellid"
metadata1$ID <- paste0("a", "_", rownames(metadata1))
rownames(metadata1) <- metadata1$ID

# Extract metadata information from Sten's DG
metadata2 <- data.frame(seurat2@meta.data)
metadata2 <- metadata2[, -19]
Cellid <- rownames(metadata2)
metadata2 <- cbind(metadata2, new_col = Cellid)
colnames(metadata2)[26] <- "Cellid"
colnames(metadata2)[23] <- "Cluster"
metadata2$ID <- paste0("b", "_", rownames(metadata2))
rownames(metadata2) <- metadata2$ID

# Merge two datasets
seurat3 <- merge(seurat1, y = seurat2, 
                 add.cell.ids = c("a", "b"), 
                 project = "DG")
######################################################################################################
# Adjust metadata
######################################################################################################
# GEMgroup
seurat3$GEMgroup <- recode(seurat3$GEMgroup, # Semisupervised annotation
                           "a1" = "1", 
                           "a2" = "2", 
                           "a3" = "3", 
                           "a4" = "4", 
                           "b1" = "5", 
                           "b2" = "6", 
                           "b3" = "7", 
                           "b4" = "8", 
                           "b5" = "9", 
                           "b6" = "10", 
                           "b7" = "11", 
                           "b8" = "12", 
                           "b9" = "13", 
                           "bb10" = "14", 
                           "b1b1" = "15")
# Batch
seurat3$Batch <- seurat3$GEMgroup
seurat3$Batch <- recode(seurat3$Batch, # Semisupervised annotation
                        "1" = "1", 
                        "2" = "1", 
                        "3" = "1", 
                        "4" = "1", 
                        "5" = "2", 
                        "6" = "2", 
                        "7" = "3", 
                        "8" = "3", 
                        "9" = "4", 
                        "10" = "5", 
                        "11" = "5", 
                        "12" = "5", 
                        "13" = "6", 
                        "14" = "6", 
                        "15" = "6")
# Dataset
seurat3$Dataset <- seurat3$GEMgroup
seurat3$Dataset <- recode(seurat3$Dataset, # Semisupervised annotation
                          "1" = "Current", 
                          "2" = "Current", 
                          "3" = "Current", 
                          "4" = "Current", 
                          "5" = "Hochgerner 2018", 
                          "6" = "Hochgerner 2018", 
                          "7" = "Hochgerner 2018", 
                          "8" = "Hochgerner 2018", 
                          "9" = "Hochgerner 2018", 
                          "10" = "Hochgerner 2018", 
                          "11" = "Hochgerner 2018", 
                          "12" = "Hochgerner 2018", 
                          "13" = "Hochgerner 2018", 
                          "14" = "Hochgerner 2018", 
                          "15" = "Hochgerner 2018")
my_levels <- c("Current", "Hochgerner 2018")
seurat3$Dataset <- factor(x = seurat3$Dataset, levels = my_levels)
######################################################################################################
# harmony
######################################################################################################
library(harmony)
seurat3 <- NormalizeData(seurat3, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat3 <- RunHarmony(seurat3, group.by.vars = "Batch", dims.use = 1:30, max.iter.harmony = 50) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
DimPlot(seurat3, reduction = "umap", label = TRUE)
DimPlot(seurat3, reduction = "tsne", label = TRUE)
# Save the metadata
metadata3 <- data.frame(seurat3@meta.data)
metadata3 <- metadata3[, c(20, 33, 34)]
######################################################################################################
# Transfer labels
######################################################################################################
metadata1 <- metadata1[rownames(metadata1) %in% colnames(seurat3), ]
metadata2 <- metadata2[rownames(metadata2) %in% colnames(seurat3), ]
seurat3 <- AddMetaData(seurat3, metadata = metadata1)
seurat3 <- AddMetaData(seurat3, metadata = metadata2)
seurat3 <- AddMetaData(seurat3, metadata = metadata3)
######################################################################################################
# Adjust labels
######################################################################################################
# Create "IntergrativeAnnotation"
## Celltype
Idents(seurat3) <- "Celltype"
Astrocyte <- subset(seurat3, idents = c("Astrocyte"))
qNSC <- subset(seurat3, idents = c("qNSC"))
nIPC_1 <- subset(seurat3, idents = c("nIPC"))
Neuroblast_1 <- subset(seurat3, idents = c("Neuroblast"))
GC_1 <- subset(seurat3, idents = c("GC"))
Cajal_Retzius_1 <- subset(seurat3, idents = c("Cajal-Retzius"))
OPC_1 <- subset(seurat3, idents = c("OPC"))
Endothelial_1 <- subset(seurat3, idents = c("Endothelial"))
Pericyte <- subset(seurat3, idents = c("Pericyte"))
SMC <- subset(seurat3, idents = c("SMC"))
Microglia_1 <- subset(seurat3, idents = c("Microglia"))
## Cluster
Idents(seurat3) <- "Cluster"
Astro <- subset(seurat3, idents = c("Astro"))
RGL <- subset(seurat3, idents = c("RGL"))
nIPC_2 <- subset(seurat3, idents = c("nIPC"))
Neuroblast_2 <- subset(seurat3, idents = c("Neuroblast"))
Immature_GC <- subset(seurat3, idents = c("Immature-GC"))
GC_2 <- subset(seurat3, idents = c("GC"))
Immature_Pyr <- subset(seurat3, idents = c("Immature-Pyr"))
CA3_Pyr <- subset(seurat3, idents = c("CA3-Pyr"))
Immature_GABA <- subset(seurat3, idents = c("Immature-GABA"))
GABA <- subset(seurat3, idents = c("GABA"))
Cajal_Retzius_2 <- subset(seurat3, idents = c("Cajal-Retzius"))
Ependymal <- subset(seurat3, idents = c("Ependymal"))
OPC_2 <- subset(seurat3, idents = c("OPC"))
NFOL <- subset(seurat3, idents = c("NFOL"))
MOL <- subset(seurat3, idents = c("MOL"))
Endothelial_2 <- subset(seurat3, idents = c("Endothelial"))
VLMC <- subset(seurat3, idents = c("VLMC"))
Microglia_2 <- subset(seurat3, idents = c("Microglia"))
PVM <- subset(seurat3, idents = c("PVM"))

Astrocyte <- WhichCells(Astrocyte)
Astro <- WhichCells(Astro)
qNSC <- WhichCells(qNSC)
RGL <- WhichCells(RGL)
nIPC_1 <- WhichCells(nIPC_1)
nIPC_2 <- WhichCells(nIPC_2)
Neuroblast_1 <- WhichCells(Neuroblast_1)
Neuroblast_2 <- WhichCells(Neuroblast_2)
Immature_GC <- WhichCells(Immature_GC)
GC_1 <- WhichCells(GC_1)
GC_2 <- WhichCells(GC_2)
Immature_Pyr <- WhichCells(Immature_Pyr)
CA3_Pyr <- WhichCells(CA3_Pyr)
Immature_GABA <- WhichCells(Immature_GABA)
GABA <- WhichCells(GABA)
Cajal_Retzius_1 <- WhichCells(Cajal_Retzius_1)
Cajal_Retzius_2 <- WhichCells(Cajal_Retzius_2)
Ependymal <- WhichCells(Ependymal)
OPC_1 <- WhichCells(OPC_1)
OPC_2 <- WhichCells(OPC_2)
NFOL <- WhichCells(NFOL)
MOL <- WhichCells(MOL)
Endothelial_1 <- WhichCells(Endothelial_1)
Endothelial_2 <- WhichCells(Endothelial_2)
Pericyte <- WhichCells(Pericyte)
SMC <- WhichCells(SMC)
VLMC <- WhichCells(VLMC)
Microglia_1 <- WhichCells(Microglia_1)
Microglia_2 <- WhichCells(Microglia_2)
PVM <- WhichCells(PVM)

seurat3$IntegrativeAnnotation <- "NA"
Idents(seurat3) <- "IntegrativeAnnotation"
Idents(seurat3, cells = Astrocyte) <- "Astrocyte" 
Idents(seurat3, cells = Astro) <- "Astrocyte" 
Idents(seurat3, cells = qNSC) <- "RGL" 
Idents(seurat3, cells = RGL) <- "RGL" 
Idents(seurat3, cells = nIPC_1) <- "nIPC" 
Idents(seurat3, cells = nIPC_2) <- "nIPC" 
Idents(seurat3, cells = Neuroblast_1) <- "Neuroblast" 
Idents(seurat3, cells = Neuroblast_2) <- "Neuroblast" 
Idents(seurat3, cells = Immature_GC) <- "GC" 
Idents(seurat3, cells = GC_1) <- "GC" 
Idents(seurat3, cells = GC_2) <- "GC" 
Idents(seurat3, cells = Immature_Pyr) <- "CA3-Pyr" 
Idents(seurat3, cells = CA3_Pyr) <- "CA3-Pyr" 
Idents(seurat3, cells = Immature_GABA) <- "GABA" 
Idents(seurat3, cells = GABA) <- "GABA" 
Idents(seurat3, cells = Cajal_Retzius_1) <- "Cajal-Retzius" 
Idents(seurat3, cells = Cajal_Retzius_2) <- "Cajal-Retzius" 
Idents(seurat3, cells = Ependymal) <- "Ependymal" 
Idents(seurat3, cells = OPC_1) <- "OPC" 
Idents(seurat3, cells = OPC_2) <- "OPC" 
Idents(seurat3, cells = NFOL) <- "Oligodendrocyte" 
Idents(seurat3, cells = MOL) <- "Oligodendrocyte" 
Idents(seurat3, cells = Endothelial_1) <- "Vascular" 
Idents(seurat3, cells = Endothelial_2) <- "Vascular" 
Idents(seurat3, cells = Pericyte) <- "Vascular" 
Idents(seurat3, cells = SMC) <- "Vascular" 
Idents(seurat3, cells = VLMC) <- "VLMC" 
Idents(seurat3, cells = Microglia_1) <- "Microglia" 
Idents(seurat3, cells = Microglia_2) <- "Microglia" 
Idents(seurat3, cells = PVM) <- "PVM" 
seurat3$IntegrativeAnnotation <- seurat3@active.ident
my_levels <- c("Astrocyte", "RGL", "nIPC", "Neuroblast", "GC", "CA3-Pyr", "GABA", "Cajal-Retzius", 
               "Ependymal", "OPC", "Oligodendrocyte", "Vascular", "VLMC", "Microglia", "PVM")
seurat3$IntegrativeAnnotation <- factor(x = seurat3$IntegrativeAnnotation, levels = my_levels)
Idents(seurat3) <- "IntegrativeAnnotation"
# Adjust sample age
seurat3$SampleTime <- recode(seurat3$GEMgroup, # Semisupervised annotation
                             "1" = "Young", 
                             "2" = "Middle-age", 
                             "3" = "Old", 
                             "4" = "Old", 
                             "5" = "Young", 
                             "6" = "Juvenile", 
                             "7" = "Young", 
                             "8" = "Young", 
                             "9" = "Perinatal", 
                             "10" = "Perinatal", 
                             "11" = "Perinatal", 
                             "12" = "Juvenile", 
                             "13" = "Perinatal", 
                             "14" = "Perinatal", 
                             "15" = "Juvenile")
my_levels <- c("Perinatal", "Juvenile", "Young", "Middle-age", "Old")
seurat3$SampleTime <- factor(x = seurat3$SampleTime, levels = my_levels)
######################################################################################################
# Visualization
######################################################################################################
DimPlot(seurat3, reduction = "umap")
DimPlot(seurat3, reduction = "umap", 
        cols = c("#dd3497", "#993404", "#fe9929", "#ef3b2c", "#ec7014", "#fb6a4a", 
                 "#969696", "#a8ddb5", "#4eb3d3", "#1d91c0", "#f7fcb9", "#9e9ac8", 
                 "#88419d", "#67000d", "#54278f", "#41ae76", "#000000"))
ggsave("Integration_Annotation_Hochgerner_UMAP.png", width = 10, height = 8)
ggsave("Integration_Annotation_Hochgerner_UMAP.pdf", width = 10, height = 8)
p1 <- DimPlot(seurat3, reduction = "umap", group.by = "Celltype", label = TRUE)
p2 <- DimPlot(seurat3, reduction = "umap", group.by = "Cluster", label = TRUE)
p3 <- DimPlot(seurat3, reduction = "umap", group.by = "SampleTime", label = TRUE)
p4 <- DimPlot(seurat3, reduction = "umap", group.by = "Dataset", label = TRUE)
p5 <-  p1|p2|p3|p4
ggsave("Integration_Hochgerner_UMAP.png", plot = p5, width = 40, height = 8)
ggsave("Integration_Hochgerner_UMAP.pdf", plot = p5, width = 40, height = 8)
######################################################################################################
# Assign cell cycle genes
######################################################################################################
CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(seurat3))
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(seurat3))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(seurat3))
seurat3 <- CellCycleScoring(object = seurat3, g2m.features = g2m_genes, s.features = s_genes)
FeaturePlot(seurat3, features = c("S.Score", "G2M.Score"), reduction = "umap", 
            min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seurat3, features = c("S.Score", "G2M.Score"), 
            min.cutoff = "q10", max.cutoff = "q90", cols = c("grey90", brewer.pal(9, "YlGnBu")))

