##== Vascular ==##

# Load data
Idents(seurat_object.integrated) <- "Celltype"
seurat_vascular <- subset(seurat_object.integrated, idents = c("Endothelial", "Pericyte", "SMC"))
######################################################################################################
# Re-clustering
######################################################################################################
DefaultAssay(seurat_vascular) <- "RNA"
seurat_vascular <- NormalizeData(seurat_vascular, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 6000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_vascular <- RunUMAP(seurat_vascular, reduction = "pca", dims = 1:20) %>% 
  RunTSNE(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
Idents(seurat_vascular) <- "Celltype"
DimPlot(seurat_vascular, reduction = "umap", cols = c("#78c679", "#6baed6", "#225ea8"))
ggsave("UMAP_Vascular.png", width = 5.5, height = 4)
ggsave("UMAP_Vascular.pdf", width = 5.5, height = 4)
# Adjust the level of metadata
my_levels <- c("Early_Endothelial", "Endothelial", "Pericyte", "SMC")
seurat_vascular$Subtype <- factor(x = seurat_vascular$Subtype, levels = my_levels)
Idents(seurat_vascular) <- "Subtype"
DimPlot(seurat_vascular, reduction = "umap", cols = c("#cb181d", "#fe9929", "#1d91c0", "#78c679"))
ggsave("UMAP_Vascular_Subtype.png", width = 5.8, height = 4)
ggsave("UMAP_Vascular_Subtype.pdf", width = 5.8, height = 4)

##== Visualization ==##
######################################################################################################
# Age distribution
######################################################################################################
Idents(seurat_vascular) <- "Age"
DimPlot(seurat_vascular, reduction = "umap", cols.highlight = "#3182bd", 
        cells.highlight = WhichCells(seurat_vascular, idents = c("Young"))) + NoLegend()
ggsave("UMAP_Vascular_Young.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Young.pdf", width = 4.4, height = 4)
DimPlot(seurat_vascular, reduction = "umap", cols.highlight = "#31a354", 
        cells.highlight = WhichCells(seurat_vascular, idents = c("Middle"))) + NoLegend()
ggsave("UMAP_Vascular_Middle.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Middle.pdf", width = 4.4, height = 4)
DimPlot(seurat_vascular, reduction = "umap", cols.highlight = "#de2d26", 
        cells.highlight = WhichCells(seurat_vascular, idents = c("Old"))) + NoLegend()
ggsave("UMAP_Vascular_Old.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Old.pdf", width = 4.4, height = 4)

######################################################################################################
# Sex distribution
######################################################################################################
Idents(seurat_vascular) <- "Sex"
DimPlot(seurat_vascular, reduction = "umap", cols.highlight = "#3182bd", 
        cells.highlight = WhichCells(seurat_vascular, idents = c("Female"))) + NoLegend()
ggsave("UMAP_Vascular_Female.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Female.pdf", width = 4.4, height = 4)
DimPlot(seurat_vascular, reduction = "umap", cols.highlight = "#de2d26", 
        cells.highlight = WhichCells(seurat_vascular, idents = c("Male"))) + NoLegend()
ggsave("UMAP_Vascular_Male.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Male.pdf", width = 4.4, height = 4)


Idents(seurat_female) <- "Age"
DimPlot(seurat_female, reduction = "umap", cols.highlight = "#3182bd", 
        cells.highlight = WhichCells(seurat_female, idents = c("Young"))) + NoLegend()
ggsave("UMAP_Vascular_Young_Female.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Young_Female.pdf", width = 4.4, height = 4)
DimPlot(seurat_female, reduction = "umap", cols.highlight = "#31a354", 
        cells.highlight = WhichCells(seurat_female, idents = c("Middle"))) + NoLegend()
ggsave("UMAP_Vascular_Middle_Female.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Middle_Female.pdf", width = 4.4, height = 4)
DimPlot(seurat_female, reduction = "umap", cols.highlight = "#de2d26", 
        cells.highlight = WhichCells(seurat_female, idents = c("Old"))) + NoLegend()
ggsave("UMAP_Vascular_Old_Female.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Old_Female.pdf", width = 4.4, height = 4)


Idents(seurat_male) <- "Age"
DimPlot(seurat_male, reduction = "umap", cols.highlight = "#3182bd", 
        cells.highlight = WhichCells(seurat_male, idents = c("Young"))) + NoLegend()
ggsave("UMAP_Vascular_Young_Male.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Young_Male.pdf", width = 4.4, height = 4)
DimPlot(seurat_male, reduction = "umap", cols.highlight = "#31a354", 
        cells.highlight = WhichCells(seurat_male, idents = c("Middle"))) + NoLegend()
ggsave("UMAP_Vascular_Middle_Male.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Middle_Male.pdf", width = 4.4, height = 4)
DimPlot(seurat_male, reduction = "umap", cols.highlight = "#de2d26", 
        cells.highlight = WhichCells(seurat_male, idents = c("Old"))) + NoLegend()
ggsave("UMAP_Vascular_Old_Male.png", width = 4.4, height = 4)
ggsave("UMAP_Vascular_Old_Male.pdf", width = 4.4, height = 4)

######################################################################################################
# FeaturePlot
######################################################################################################
DefaultAssay(seurat_vascular) <- "RNA"
FeaturePlot(seurat_vascular, reduction = "umap", 
            features = c("Cldn5", "Pecam1", "Flt1", 
                         "Vtn", "Rgs5", "Pdgfrb", 
                         "Cnn1", "Tagln", "Acta2"), 
            pt.size = 0.5, cols = c("lightgrey", "#ae017e"), 
            min.cutoff = "q10", max.cutoff = "q90")
ggsave("Featureplot_Vascular.png", width = 13.5, height = 12)
ggsave("Featureplot_Vascular.pdf", width = 13.5, height = 12)

######################################################################################################
# Dotplot
######################################################################################################
Idents(seurat_vascular) <- "Celltype"
DotPlot(seurat_vascular, 
        features = c("Cldn5", "Pecam1", "Flt1", 
                     "Vtn", "Rgs5", "Pdgfrb", 
                     "Cnn1", "Tagln", "Acta2"),  
        cols = c("#f768a1", "#ae017e", "#49006a"), dot.scale = 10) + 
  RotatedAxis() + 
  coord_flip() + 
  theme(panel.background = element_rect(size = 1.5, fill = "white", colour = "black")) 
ggsave("Dotplot_Vascular.png", width = 6, height = 5)
ggsave("Dotplot_Vascular.pdf", width = 6, height = 5)




















##== Microglia ==##

# Load data
Idents(seurat_object.integrated) <- "Celltype"
seurat_immune <- subset(seurat_object.integrated, idents = c("Microglia"))

######################################################################################################
# Re-clustering
######################################################################################################
DefaultAssay(seurat_immune) <- "RNA"
seurat_immune <- NormalizeData(seurat_immune, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_immune <- RunUMAP(seurat_immune, reduction = "pca", dims = 1:15) %>% 
  RunTSNE(reduction = "pca", dims = 1:15) %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = 0.5)
Idents(seurat_immune) <- "seurat_clusters"
DimPlot(seurat_immune, reduction = "umap", label = T)

######################################################################################################
# Expression of exAM genes
######################################################################################################
# Read exAM list from Marsh 2022
exAM <- c("Fos", "Junb", "Zfp36", "Jun", "Hspa1a", "Socs3", "Rgs1", "Egr1", "Btg2", "Fosb", "Hist1h1d", 
          "Ier5", "1500015O10Rik", "Atf3", "Hist1h2ac", "Dusp1", "Hist1h1e", "Folr1", "Serpine1")

# FeaturePlot
DimPlot(seurat_immune, reduction = "umap", label = T, group.by = "seurat_clusters")
FeaturePlot(seurat_immune, 
            features = exAM, 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, 
            cols = c("lightgrey", "Darkred"), ncol = 3)
ggsave("FeaturePlot_Immune_exAM.png", width = 13.5, height = 24)
ggsave("FeaturePlot_Immune_exAM.pdf", width = 13.5, height = 24)

# Module score
DefaultAssay(seurat_immune) <- "RNA"
seurat_immune <- AddModuleScore(seurat_immune, 
                                features = list(exAM), 
                                name = "exAM", 
                                assay = "RNA")
seurat_immune$exAM <- seurat_immune$exAM1
seurat_immune$exAM1 <- NULL
FeaturePlot(seurat_immune, features = c("exAM"), 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, 
            cols = c("lightgrey", "Darkred"), ncol = 1)
ggsave("FeaturePlot_Immune_exAM_ModuleScore.png", width = 4.5, height = 4)
ggsave("FeaturePlot_Immune_exAM_ModuleScore.pdf", width = 4.5, height = 4)

ggplot(seurat_immune@meta.data, aes(x = exAM)) + 
  geom_histogram(bins = 100) + 
  geom_vline(aes(xintercept = 0), color  = "red", linetype = 2) + 
  theme_bw()
ggsave("Histogram_exAM.png", width = 4.5, height = 4.5)
ggsave("Histogram_exAM.pdf", width = 4.5, height = 4.5)

# Update SeuratObject
saveRDS(seurat_immune, "seurat_immune.rds")

######################################################################################################
# Regress out exAM genes
######################################################################################################
exprMat <- as.matrix(seurat_immune@assays$RNA@counts)
exprMat <- exprMat[!rownames(exprMat) %in% exAM, ]
dim(exprMat)
meta_immune <- as.data.frame(seurat_immune@meta.data)

# Re-clustering
seurat <- CreateSeuratObject(counts = exprMat)
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:15) %>% 
  RunTSNE(reduction = "pca", dims = 1:15) %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = 0.5)
seurat <- AddMetaData(seurat, metadata = meta_immune)
Idents(seurat) <- "Subtype"
DimPlot(seurat, reduction = "umap", cols = c("#fe9929", "#ef3b2c", "#525252"))
# Adjust the level of metadata
my_levels <- c("Microglia_1", "Microglia_2", "Microglia_3")
seurat$Subtype <- factor(x = seurat$Subtype, levels = my_levels)
Idents(seurat) <- "Subtype"
DimPlot(seurat, reduction = "umap", cols = c("#9e9ac8", "#f768a1", "#737373"))
ggsave("UMAP_Microglia_no_exAM.png", width = 5.5, height = 4)
ggsave("UMAP_Microglia_no_exAM.pdf", width = 5.5, height = 4)

######################################################################################################
# Save data
######################################################################################################
saveRDS(seurat, "seurat_immune_no_exAM.rds")







##== Visualization ==##
######################################################################################################
# Age distribution
######################################################################################################
Idents(seurat_immune) <- "Age"
DimPlot(seurat_immune, reduction = "umap", cols.highlight = "#3182bd", 
        cells.highlight = WhichCells(seurat_immune, idents = "Young"))
ggsave("UMAP_Immune_Young.png", width = 5.5, height = 4)
ggsave("UMAP_Immune_Young.pdf", width = 5.5, height = 4)
DimPlot(seurat_immune, reduction = "umap", cols.highlight = "#31a354", 
        cells.highlight = WhichCells(seurat_immune, idents = "Middle"))
ggsave("UMAP_Immune_Middle.png", width = 5.5, height = 4)
ggsave("UMAP_Immune_Middle.pdf", width = 5.5, height = 4)
DimPlot(seurat_immune, reduction = "umap", cols.highlight = "#de2d26", 
        cells.highlight = WhichCells(seurat_immune, idents = "Old"))
ggsave("UMAP_Immune_Old.png", width = 5.5, height = 4)
ggsave("UMAP_Immune_Old.pdf", width = 5.5, height = 4)







