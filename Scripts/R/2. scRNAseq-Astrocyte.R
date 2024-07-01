Idents(seurat_object.integrated) <- "Celltype"
seurat_AS <- subset(seurat_object.integrated, idents = c("Astrocyte"))

######################################################################################################
# Re-clustering
######################################################################################################
DefaultAssay(seurat_AS) <- "RNA"
seurat_AS <- NormalizeData(seurat_AS, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_AS <- RunUMAP(seurat_AS, reduction = "pca", dims = 1:30) %>% 
  RunTSNE(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.5)
DimPlot(seurat_AS, reduction = "umap", label = T)

######################################################################################################
# Annotation
######################################################################################################
seurat_AS$Subtype <- recode(seurat_AS$seurat_clusters, 
                            "0" = "Astrocyte_1", 
                            "1" = "Astrocyte_2", 
                            "2" = "Astrocyte_1", 
                            "3" = "Astrocyte_1", 
                            "4" = "Astrocyte_2", 
                            "5" = "Astrocyte_1", 
                            "6" = "Astrocyte_1")
# Adjust the level of metadata
my_levels <- c("Astrocyte_1", "Astrocyte_2")
seurat_AS$Subtype <- factor(x = seurat_AS$Subtype, levels = my_levels)
Idents(seurat_AS) <- "Subtype"
DimPlot(seurat_AS, reduction = "umap", cols = c("#ae017e", "#fa9fb5"))
ggsave("UMAP_Astrocyte_New.png", width = 5.5, height = 4)
ggsave("UMAP_Astrocyte_New.pdf", width = 5.5, height = 4)

######################################################################################################
# Dotplot
######################################################################################################
Idents(seurat_AS) <- "Subtype"
DotPlot(seurat_AS, 
        features = c("Aldh1a1", "Thrsp", "Kcnk1", 
                     "Sparc", "Nnat", "Il33"),  
        cols = c("#f768a1", "#ae017e", "#49006a"), dot.scale = 10) + 
  RotatedAxis() + 
  coord_flip() + 
  theme(panel.background = element_rect(size = 1.5, fill = "white", colour = "black")) 
ggsave("Dotplot_Astrocyte.png", width = 4.8, height = 5)
ggsave("Dotplot_Astrocyte.pdf", width = 4.8, height = 5)

######################################################################################################
# FeaturePlot
######################################################################################################
FeaturePlot(seurat_AS, features = c("Aldh1a1", "Thrsp", "Kcnk1", 
                                    "Sparc", "Nnat", "Il33"), 
            reduction = "umap", cols = c("lightgrey", "#ae017e"), 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 2, 
            ncol = 3)
ggsave("Featureplot_AS.png", width = 13.5, height = 8)
ggsave("Featureplot_AS.pdf", width = 13.5, height = 8)





######################################################################################################
# Differential expression 
######################################################################################################
# Adjust the identity
Idents(seurat_AS) <- seurat_AS$Subtype
# Differential expression
All_AS_markers <- FindMarkers(object = seurat_AS, 
                              assay = "RNA", slot = "data", 
                              ident.1 = "Astrocyte_1", 
                              ident.2 = "Astrocyte_2", 
                              only.pos = FALSE, 
                              min.pct = 0, 
                              logfc.threshold = 0, 
                              test.use = "MAST")
write.csv(All_AS_markers, "All_AS_markers.csv")
# Subtype-specific genes
Astrocyte_1 <- data %>% 
  filter(Change == "Up")
Genelist_Astrocyte_1 <- rownames(Astrocyte_1)
Astrocyte_2 <- data %>% 
  filter(Change == "Down")
Genelist_Astrocyte_2 <- rownames(Astrocyte_2)
write.csv(Genelist_Astrocyte_1, "Genelist_Astrocyte_1.csv") ## Used for GSEA
write.csv(Genelist_Astrocyte_2, "Genelist_Astrocyte_2.csv") ## Used for GSEA

# Violin plot
VlnPlot(seurat_AS, features = c("Slc1a3", "Sox2", "S100b", 
                                "Aldh1a1", "Thrsp", "Kcnk1", 
                                "Sparc", "Nnat", "Il33"), 
        pt.size = 0, ncol = 3, cols = c("#ae017e", "#fa9fb5")) 
ggsave("Vlnplot_AS_Subtype.png", width = 7.5, height = 12)
ggsave("Vlnplot_AS_Subtype.pdf", width = 7.5, height = 12)

