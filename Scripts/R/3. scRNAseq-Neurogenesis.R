##== Whole neurogenic lineage ==##

Idents(seurat_object.integrated) <- "Celltype"
seurat_neurogenic <- subset(seurat_object.integrated, idents = c("qNSC", "nIPC", "Neuroblast"))
######################################################################################################
# Re-clustering
######################################################################################################
DefaultAssay(seurat_neurogenic) <- "RNA"
seurat_neurogenic <- NormalizeData(seurat_neurogenic, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_neurogenic <- RunUMAP(seurat_neurogenic, reduction = "pca", dims = 1:12) %>% 
  RunTSNE(reduction = "pca", dims = 1:12) %>% 
  FindNeighbors(reduction = "pca", dims = 1:12) %>% 
  FindClusters(resolution = 1.0)
# Adjust the level of metadata
my_levels <- c("qNSC", "nIPC", "Neuroblast")
seurat_neurogenic$Celltype <- factor(x = seurat_neurogenic$Celltype, levels = my_levels)
Idents(seurat_neurogenic) <- "Celltype"
DimPlot(seurat_neurogenic, reduction = "umap", cols = c("#993404", "#fe9929", "#ef3b2c"), 
        pt.size = 2)
ggsave("UMAP_Neurogenic.png", width = 5.5, height = 4)
ggsave("UMAP_Neurogenic.pdf", width = 5.5, height = 4)

######################################################################################################
# Feature gene expression
######################################################################################################
# Dotplot
selected.genes <- c("Gfap", "Rfx4", "Ascl1", 
                    "Top2a", "Neurog2", 
                    "Neurod1", "Neurod2", "Calb2")
DotPlot(seurat_neurogenic, 
        features = selected.genes, 
        cols = c("#f768a1", "#ae017e", "#49006a"), dot.scale = 10) + 
  RotatedAxis() + 
  coord_flip() + 
  theme(panel.background = element_rect(size = 1.5, fill = "white", colour = "black")) 
ggsave("Dotplot_Neurogenic.png", width = 6, height = 5)
ggsave("Dotplot_Neurogenic.pdf", width = 6, height = 5)

# Featureplot
FeaturePlot(seurat_neurogenic, features = c("Gfap", "Ascl1", "Top2a", "Neurog2", "Neurod1", "Calb2"), 
            reduction = "umap", cols = c("lightgrey", "#ae017e"), 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 2, 
            ncol = 3)
ggsave("Featureplot_Neurogenic.png", width = 13.5, height = 8)
ggsave("Featureplot_Neurogenic.pdf", width = 13.5, height = 8)

FeaturePlot(seurat_neurogenic, features = c("S.Score", "G2M.Score"), 
            reduction = "umap", cols = c("lightgrey", "#ae017e"), 
            min.cutoff = "q10", max.cutoff = "q90", pt.size = 2, 
            ncol = 1)
ggsave("CellcycleScore_neurogenic.png", width = 4.5, height = 8)
ggsave("CellcycleScore_neurogenic.pdf", width = 4.5, height = 8)










##== SCENIC ==##
##== On server ==##
##== Data preparation ==##
######################################################################################################
# Load library 
######################################################################################################
rm(list=ls())
options(stringsAsFactors = FALSE) # Do not convert strings to factors
library(Seurat)
library(tidyverse)
library(patchwork)
library(parallel)
library(AUCell)
library(GENIE3)
library(RcisTarget)
library(SCENIC)

######################################################################################################
# Load data
######################################################################################################
scRNAseq <- readRDS("./Data/SeuratObject/Total/Final_Annotation.rds")
scRNAseq <- SetIdent(scRNAseq, value = "Celltype")
DefaultAssay(scRNAseq) <- "RNA"

# Set up working path (fit the path on server)
dir.create("./SCENIC")
dir.create("./SCENIC/int")
setwd("./SCENIC") 

# Prepare metadata
cellInfo <- data.frame(scRNAseq@meta.data)
colnames(cellInfo)[which(colnames(cellInfo) == "Age")] <- "Age"
colnames(cellInfo)[which(colnames(cellInfo) == "Celltype")] <- "Celltype"
cellInfo <- cellInfo[, c("Age", "Celltype")]
saveRDS(cellInfo, file = "./int/cellInfo.Rds")

# Prepare expression matrix
exprMat <- as.matrix(scRNAseq@assays$RNA@data) # No big difference between raw and normalized counts

######################################################################################################
# Set up the environment for calculation
######################################################################################################
mydbDIR <- "~/Desktop/scRNAseq/Projects/Common/cisTarget_databases/"
mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather", 
           "mm9-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org = "mgi", 
                                  nCores = detectCores(),
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "Neurogenesis")
saveRDS(scenicOptions, "int/scenicOptions.rds")
saveRDS(exprMat, "int/exprMat.rds") # For the use on local

genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions, nParts = 3)


##== On local ==##
# Set up working path
setwd("./SCENIC") 
exprMat <- readRDS("int/exprMat.rds")
scenicOptions <- readRDS("int/scenicOptions.rds")
scenicOptions@settings$dbDir <- "~/Desktop/scRNAseq/Projects/Common/cisTarget_databases/" # RcisTarget databases location
scenicOptions@settings$nCores <- 4

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_binary_nonDupl") 
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat)
savedSelections <- shiny::runApp(aucellApp)

newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file = getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

runSCENIC_4_aucell_binarize(scenicOptions, exprMat)











##== monocle2 ==##
######################################################################################################
# Getting Started with Monocle
######################################################################################################
# Importing data from Seurat object
data <- as(as.matrix(seurat_neurogenic@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = seurat_neurogenic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

# Create CellDataSet object
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, 
                             relative_expr = T, 
                             cores = detectCores())

# Filtering low-quality cells
mycds <- detectGenes(mycds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(mycds), 
                                    num_cells_expressed >= 10)) 

######################################################################################################
# Constructing single cell trajectories
######################################################################################################
# Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(mycds[expressed_genes, ], 
                                      fullModelFormulaStr = "~Celltype", 
                                      cores = detectCores()) 
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) # 6205 genes

mycds <- setOrderingFilter(mycds, ordering_genes)  
plot_ordering_genes(mycds)

# Trajectory step 2: reduce data dimensionality
mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree")

# Trajectory step 3: order cells along the trajectory
mycds <- orderCells(mycds, reverse = T)

######################################################################################################
# Differential Expression Analysis
######################################################################################################
# Clustering Genes by Pseudotemporal Expression Pattern
diff_test_res <- differentialGeneTest(mycds[ordering_genes, ],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01)) # 5415 genes
p <- plot_pseudotime_heatmap(mycds[sig_gene_names, ],
                             num_clusters = 4,
                             cores = detectCores(),
                             show_rownames = T, 
                             hmcols = viridis(1000, option = "inferno", alpha = 1, begin = 0), 
                             return_heatmap = T)
ggsave("Monocle_Heatmap_Neurogenic.png", p, width = 8, height = 10)
ggsave("Monocle_Heatmap_Neurogenic.pdf", p, width = 8, height = 10)

# Extract pseudotemporally ordered genes in qNSCs
t <- as.data.frame(cutree(p$tree_row, k = 4))
colnames(t) <- "Cluster"
t$Cluster <- recode(t$Cluster, 
                    "1" = "module_1", 
                    "4" = "module_2", 
                    "3" = "module_3", 
                    "2" = "module_4")
t$Gene <- rownames(t)
module_1 <- t %>% 
  filter(Cluster == "module_1")
module_1 <- rownames(module_1)
module_2 <- t %>% 
  filter(Cluster == "module_2")
module_2 <- rownames(module_2)
module_3 <- t %>% 
  filter(Cluster == "module_3")
module_3 <- rownames(module_3)
module_4 <- t %>% 
  filter(Cluster == "module_4")
module_4 <- rownames(module_4)
write.csv(t, "pseudotime.Neurogenic.All.csv") # Used for GSEA
write.csv(module_1, "pseudotime.Neurogenic.module.1.csv") # Used for GSEA
write.csv(module_2, "pseudotime.Neurogenic.module.2.csv") # Used for GSEA
write.csv(module_3, "pseudotime.Neurogenic.module.3.csv") # Used for GSEA
write.csv(module_4, "pseudotime.Neurogenic.module.4.csv") # Used for GSEA

# Asign module score
DefaultAssay(seurat_neurogenic) <- "RNA"
seurat_neurogenic <- AddModuleScore(seurat_neurogenic, 
                                    features = list(module_1), 
                                    name = "module_1", 
                                    assay = "RNA")
seurat_neurogenic <- AddModuleScore(seurat_neurogenic, 
                                    features = list(module_2), 
                                    name = "module_2", 
                                    assay = "RNA")
seurat_neurogenic <- AddModuleScore(seurat_neurogenic, 
                                    features = list(module_3), 
                                    name = "module_3", 
                                    assay = "RNA")
seurat_neurogenic <- AddModuleScore(seurat_neurogenic, 
                                    features = list(module_4), 
                                    name = "module_4", 
                                    assay = "RNA")
seurat_neurogenic$module_1 <- seurat_neurogenic$module_11
seurat_neurogenic$module_2 <- seurat_neurogenic$module_21
seurat_neurogenic$module_3 <- seurat_neurogenic$module_31
seurat_neurogenic$module_4 <- seurat_neurogenic$module_41
seurat_neurogenic$module_11 <- NULL
seurat_neurogenic$module_21 <- NULL
seurat_neurogenic$module_31 <- NULL
seurat_neurogenic$module_41 <- NULL
Idents(seurat_neurogenic) <- "Celltype"
VlnPlot(seurat_neurogenic, features = c("module_1", "module_2", "module_3", "module_4"), 
        cols = c("#993404", "#fe9929", "#3690c0"), pt.size = 0, ncol = 4)
ggsave("Vlnplot_Pseudotemporal_gene_Neurogenic.png", width = 13.5, height = 4.5)
ggsave("Vlnplot_Pseudotemporal_gene_Neurogenic.pdf", width = 13.5, height = 4.5)
FeaturePlot(seurat_neurogenic, features = c("module_1", "module_2", "module_3", "module_4"), 
            reduction = "monocle",min.cutoff = "q10", max.cutoff = "q90", pt.size = 1.2, 
            cols = c(rgb(0.7, 0.7, 0.7, alpha = 0.2), "darkred"), ncol = 4)
ggsave("Monocle_Featureplot_Module_Neurogenic.png", width = 18, height = 4)
ggsave("Monocle_Featureplot_Module_Neurogenic.pdf", width = 18, height = 4)






##== qNSC ==##

Idents(seurat_neurogenic) <- "Celltype"
seurat_qNSC <- subset(seurat_neurogenic, ident = "qNSC")
##== Data process ==##
######################################################################################################
# Getting Started with Monocle
######################################################################################################
# Importing data from Seurat object
data <- as(as.matrix(seurat_qNSC@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = seurat_qNSC@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

# Create CellDataSet object
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, 
                             relative_expr = T, 
                             cores = detectCores())

# Filtering low-quality cells
mycds <- detectGenes(mycds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(mycds), 
                                    num_cells_expressed >= 10)) 

######################################################################################################
# Constructing single cell trajectories
######################################################################################################
# Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(mycds[expressed_genes, ], 
                                      fullModelFormulaStr = "~Age", 
                                      cores = detectCores()) 
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) # 118 genes

mycds <- setOrderingFilter(mycds, ordering_genes)  
plot_ordering_genes(mycds)

# Trajectory step 2: reduce data dimensionality
mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree")

# Trajectory step 3: order cells along the trajectory
mycds <- orderCells(mycds, reverse = T)

######################################################################################################
# Differential Expression Analysis
######################################################################################################
# Clustering Genes by Pseudotemporal Expression Pattern
diff_test_res <- differentialGeneTest(mycds[ordering_genes, ],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01)) # 99 genes
p <- plot_pseudotime_heatmap(mycds[sig_gene_names, ],
                             num_clusters = 3,
                             cores = detectCores(),
                             show_rownames = T, 
                             hmcols = viridis(1000, option = "inferno", alpha = 1, begin = 0), 
                             return_heatmap = T)
ggsave("Monocle_Heatmap_qNSC.png", p, width = 8, height = 10)
ggsave("Monocle_Heatmap_qNSC.pdf", p, width = 8, height = 10)

# Extract pseudotemporally ordered genes in qNSCs
t <- as.data.frame(cutree(p$tree_row, k = 3))
saveRDS(t, "Pseudotemporal_genes_qNSC.rds")
colnames(t) <- "Cluster"
t$Gene <- rownames(t)
module_1 <- t %>% 
  filter(Cluster == 1)
module_1 <- rownames(module_1)
module_2 <- t %>% 
  filter(Cluster == 2)
module_2 <- rownames(module_2)
module_3 <- t %>% 
  filter(Cluster == 3)
module_3 <- rownames(module_3)
write.csv(module_1, "pseudotime.qNSC.module.1.csv") # Used for GSEA
write.csv(module_2, "pseudotime.qNSC.module.2.csv") # Used for GSEA
write.csv(module_3, "pseudotime.qNSC.module.3.csv") # Used for GSEA

# Asign module score
DefaultAssay(seurat_qNSC) <- "RNA"
seurat_qNSC <- AddModuleScore(seurat_qNSC, 
                              features = list(module_1), 
                              name = "module_1", 
                              assay = "RNA")
seurat_qNSC <- AddModuleScore(seurat_qNSC, 
                              features = list(module_2), 
                              name = "module_2", 
                              assay = "RNA")
seurat_qNSC <- AddModuleScore(seurat_qNSC, 
                              features = list(module_3), 
                              name = "module_3", 
                              assay = "RNA")
seurat_qNSC$module_1 <- seurat_qNSC$module_11
seurat_qNSC$module_2 <- seurat_qNSC$module_21
seurat_qNSC$module_3 <- seurat_qNSC$module_31
seurat_qNSC$module_11 <- NULL
seurat_qNSC$module_21 <- NULL
seurat_qNSC$module_31 <- NULL
Idents(seurat_qNSC) <- "Age"
VlnPlot(seurat_qNSC, features = c("module_1", "module_2", "module_3"), 
        cols = c("#deebf7", "#6baed6", "#08519c"), pt.size = 0)
ggsave("Vlnplot_Pseudotemporal_gene_qNSC.png", width = 10, height = 4.5)
ggsave("Vlnplot_Pseudotemporal_gene_qNSC.pdf", width = 10, height = 4.5)
FeaturePlot(seurat_qNSC, features = c("module_1", "module_2", "module_3"), 
            reduction = "monocle",min.cutoff = "q10", max.cutoff = "q90", pt.size = 1.2, 
            cols = c(rgb(0.7, 0.7, 0.7, alpha = 0.2), "darkred"), ncol = 3)
ggsave("Monocle_Featureplot_Module_Neurogenic.png", width = 18, height = 4)
ggsave("Monocle_Featureplot_Module_Neurogenic.pdf", width = 18, height = 4)




##== Mfuzz ==##
######################################################################################################
# Average expression
######################################################################################################
Idents(seurat_neurogenic) <- "Celltype"
seurat_qNSC <- subset(seurat_neurogenic, idents = "qNSC")
VariableFeatures <- seurat_neurogenic@assays[["RNA"]]@var.features
avg.seurat_qNSC <- as.data.frame(log1p(AverageExpression(seurat_qNSC, group.by = "Age", verbose = F)$RNA))
avg.seurat_qNSC <- avg.seurat_qNSC[rownames(avg.seurat_qNSC) %in% VariableFeatures, ]
######################################################################################################
# Correlation
######################################################################################################
Cor_gene <- names(tail(sort(apply(avg.seurat_qNSC, 1, sd)), 1000))
pheatmap(cor(avg.seurat_qNSC[Cor_gene, ], method = "pearson"), 
         color = colorRampPalette(brewer.pal(9, "Blues"))(100), 
         show_colnames = F, border_color = NA)
######################################################################################################
# Mfuzz
######################################################################################################
dat <- as.matrix(avg.seurat_qNSC)
dat <- new("ExpressionSet", exprs = dat)

# Filtering missing values
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = "mean")
# Filtering genes which are expressed at low levels or show only small changes in expression
dat <- filter.std(dat, min.std  = 0) # 17233 genes excluded

# Standarlization
dat <- standardise(dat)

# Soft clustering of gene expression data
n <- 8
m <- mestimate(dat)
set.seed(1000)
cl <- mfuzz(dat, c = n, m = m)

# Plotting
## par("mar")
## par(mar = c(1, 1, 1, 1))

color.2 <- colorRampPalette(c("#e7e1ef", "#e7298a", "#67001f"))(100)
mfuzz.plot2(dat, cl = cl, colo = color.2, 
            mfrow = c(4, 2), centre = TRUE, x11 = F, centre.lwd = 0.2)
pdf("Mfuzz_qNSC_HVG.pdf", width = 6, height = 8)
mfuzz.plot2(dat, cl = cl, colo = color.2, 
            mfrow = c(4, 2), centre = TRUE, x11 = F, centre.lwd = 0.2)
dev.off() 
######################################################################################################
# Extract cluster information
######################################################################################################
cl$size
## [1] 666 554 519 578 451 786 438 553
gene_cluster <- cbind(cl$cluster, cl$membership)
colnames(gene_cluster)[1] <- "cluster"
gene_cluster <- as.data.frame(gene_cluster)
gene_cluster$gene <- rownames(gene_cluster)
write.csv(gene_cluster, "Mfuzz_qNSC.csv") 

######################################################################################################
# Intersection between Mfuzz and monocle 
######################################################################################################
# qNSC
gene_cluster_monocle <- read.csv("Data/DEG/Chromium/Monocle2/GeneList/pseudotime.Neurogenic.All.csv")
gene_cluster_monocle <- gene_cluster_monocle[, -1]
colnames(gene_cluster_monocle)[1] <- "cluster"
colnames(gene_cluster_monocle)[2] <- "gene"
gene_cluster_monocle$cluster <- paste0("Neurogenic", "_", gene_cluster_monocle$cluster)
gene_cluster_monocle$cluster <- recode(gene_cluster_monocle$cluster, 
                                       "Neurogenic_module_1" = "Neurogenic_1", 
                                       "Neurogenic_module_2" = "Neurogenic_2", 
                                       "Neurogenic_module_3" = "Neurogenic_3", 
                                       "Neurogenic_module_4" = "Neurogenic_4")

gene_cluster_mfuzz <- read.csv("Data/DEG/Mfuzz/GeneList/Mfuzz_qNSC.csv")
gene_cluster_mfuzz <- gene_cluster_mfuzz[, c(2, 11)]
gene_cluster_mfuzz$cluster <- paste0("qNSC", "_", gene_cluster_mfuzz$cluster)

gene_cluster_mfuzz <- gene_cluster_mfuzz %>% 
  filter(cluster != "qNSC_1") %>% 
  filter(cluster != "qNSC_3") %>% 
  filter(cluster != "qNSC_7")

gene_merge <- merge(gene_cluster_monocle, gene_cluster_mfuzz, by = "gene")
colnames(gene_merge)[2] <- "Monocle_cluster"
colnames(gene_merge)[3] <- "Mfuzz_cluster"
gene_merge <- gene_merge[, -1]

gene_merge %>% table() -> tit2d
tit <- melt(tit2d, id.vars = c("clusters", "V3"))

# Level
Neurogenic_level <- c("Neurogenic_1", "Neurogenic_2", "Neurogenic_3", "Neurogenic_4")
Neurogenic_level <- rev(Neurogenic_level)
qNSC_level <- c("qNSC_1", "qNSC_2", "qNSC_4", "qNSC_5", "qNSC_6", "qNSC_8")
qNSC_level <- rev(qNSC_level)
tit$Monocle_cluster <- ordered(tit$Monocle_cluster, levels = Neurogenic_level)
tit$Mfuzz_cluster <- ordered(tit$Mfuzz_cluster, levels = qNSC_level)
# Plotting
## Color palette
Neurogenic_1_col <- "#662506"
Neurogenic_2_col <- "#cc4c02"
Neurogenic_3_col <- "#fec44f"
Neurogenic_4_col <- "#ef3b2c"

alluvial(tit[, 2:1], 
         freq = tit$value,
         col = c(Neurogenic_1_col, Neurogenic_2_col, Neurogenic_3_col, Neurogenic_4_col), 
         border = c(Neurogenic_1_col, Neurogenic_2_col, Neurogenic_3_col, Neurogenic_4_col), 
         blocks = T, 
         hide = tit$value < 10, 
         alpha = 0.75, 
         cex = 0.7)
## Save plot
pdf("SankyPlot_qNSC_Monocle_Mfuzz.pdf")
alluvial(tit[, 2:1], 
         freq = tit$value,
         col = c(Neurogenic_1_col, Neurogenic_2_col, Neurogenic_3_col, Neurogenic_4_col), 
         border = c(Neurogenic_1_col, Neurogenic_2_col, Neurogenic_3_col, Neurogenic_4_col), 
         blocks = T, 
         hide = tit$value < 10, 
         alpha = 0.75, 
         cex = 0.7)
dev.off()

######################################################################################################
# Neurogenic Aging Signature
######################################################################################################
# Prepare genelist
Neurogenic_1 <- gene_cluster_monocle %>% 
  filter(cluster == "Neurogenic_1")
Neurogenic_3 <- gene_cluster_monocle %>% 
  filter(cluster == "Neurogenic_3")
Neurogenic_4 <- gene_cluster_monocle %>% 
  filter(cluster == "Neurogenic_4")
## Merge the age-upregulation clusters
qNSC_2 <- gene_cluster_mfuzz %>% 
  filter(cluster == "qNSC_2")
qNSC_2 <- qNSC_2[, 2]
qNSC_4 <- gene_cluster_mfuzz %>% 
  filter(cluster == "qNSC_4")
qNSC_4 <- qNSC_4[, 2]
qNSC_Up <- c(qNSC_2, qNSC_4) %>% 
  make.unique()
qNSC_Up <- qNSC_Up[qNSC_Up %in% Neurogenic_1$gene] # 354 genes
## Merge the age-downregulation clusters
qNSC_5 <- gene_cluster_mfuzz %>% 
  filter(cluster == "qNSC_5")
qNSC_5 <- qNSC_5[, 2]
qNSC_6 <- gene_cluster_mfuzz %>% 
  filter(cluster == "qNSC_6")
qNSC_6 <- qNSC_6[, 2]
qNSC_8 <- gene_cluster_mfuzz %>% 
  filter(cluster == "qNSC_8")
qNSC_8 <- qNSC_8[, 2]
qNSC_Down <- c(qNSC_5, qNSC_6, qNSC_8) %>% 
  make.unique()
qNSC_Down_1 <- qNSC_Down[qNSC_Down %in% Neurogenic_3$gene] # 354 genes
qNSC_Down_2 <- qNSC_Down[qNSC_Down %in% Neurogenic_4$gene] # 293 genes
qNSC_Down <- c(qNSC_Down_1, qNSC_Down_2) # 647 genes
# Save data
write.csv(qNSC_Up, "Mfuzz_Monocle_qNSC_Up.csv") 
write.csv(qNSC_Down, "Mfuzz_Monocle_qNSC_Down.csv") 

##== Regression ==##
library(caret)
######################################################################################################
# Get PCs
######################################################################################################
# Select pcs by % variance explained
eigs <- (seurat_neurogenic@reductions$pca@stdev)^2
PercentVAr <- data.frame(PercentVAr = eigs*100 / sum(eigs),
                         PCs = 1:length(eigs))

ggplot(PercentVAr, aes(x = as.factor(PCs), y = PercentVAr)) + 
  geom_point() + 
  geom_hline(yintercept = 2, linetype = 2, colour = "red") + 
  theme_classic()
ggsave("Regression_Model_PercentVAr.png", width = 6, height = 4)
ggsave("Regression_Model_PercentVAr.pdf", width = 6, height = 4)

# Keep only PCs explaining > 2% of the variance over the 7 first
PCs <- PercentVAr %>% filter(PercentVAr > 2) %>% pull(PCs)

# Fit the principal curve over the first 7 PCs
data <- as.data.frame(seurat_neurogenic@reductions$pca@cell.embeddings[, PCs])

fit <- principal_curve(as.matrix(data[, PCs]),
                       smoother = 'lowess',
                       trace = TRUE,
                       f = 0.5,
                       stretch = 0)
PseudoDifferentiation <- fit$lambda/max(fit$lambda)
pc.line <- as.data.frame(fit$s[order(fit$lambda), ])

data$Cluster <- seurat_neurogenic$Celltype
data$PseudoDifferentiation <- PseudoDifferentiation
## seurat_neurogenic@meta.data$PseudoDifferentiation.score <- PseudoDifferentiation.score

ggplot(data, aes(PC_1, PC_2)) + 
  geom_point(aes(color = Cluster), size = 3, shape = 16) +
  geom_line(data = pc.line, color = 'red', size = 0.77) +
  scale_color_manual(values = c("#bd0026", "#9e9ac8", "#2171b5")) + 
  theme_classic()

# Plot Pseudo-differentiation score
cols <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(100)
ggplot(data, aes(PC_1, PC_2)) +
  geom_point(aes(color = PseudoDifferentiation), size = 3, shape = 16) + 
  scale_color_gradientn(colours = rev(cols), name = 'Pseudotime score') +
  geom_line(data = pc.line, color = 'red', size = 0.77) + 
  theme_classic()

ggplot(data, aes(x = Cluster, y = PseudoDifferentiation, fill = Cluster)) +
  geom_boxplot() +
  scale_color_manual(aesthetics = "fill", values = c("#993404", "#fe9929", "#ef3b2c")) + 
  ## geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  ggtitle("Differentiation model") +
  xlab("") + 
  ylab("Differentiation score") + 
  theme_classic() + 
  NoLegend()
ggsave("Pseudo_differentiation_Neurogenic.png", width = 4, height = 4)
ggsave("Pseudo_differentiation_Neurogenic.pdf", width = 4, height = 4)

######################################################################################################
# Train a random forest regression model on the entire neurogenic lineage
######################################################################################################
DefaultAssay(seurat_qNSC) <- "RNA"
seurat_qNSC <- FindVariableFeatures(seurat_qNSC, nfeatures = 5000)

# We trained the model using only genes found to be highly variable in qNSC
model.data <- data.frame(PseudoDifferentiation = data$PseudoDifferentiation)
genes <- t(seurat_neurogenic@assays[["RNA"]]@scale.data[rownames(seurat_neurogenic@assays[["RNA"]]@scale.data) %in% 
                                                          seurat_qNSC@assays[["RNA"]]@var.features, ])
genes <- as.data.frame(as.matrix(genes))
model.data <- cbind(model.data, genes)

# We did not perform cross-validatin but relied on the Out Of Bag (OOB) error
trctrl <- trainControl(method = "none")

# Train the random forest model
set.seed(1234)
Differentiation.model <- train(PseudoDifferentiation~., 
                               data = model.data, 
                               method = "ranger",
                               trControl = trctrl,
                               importance = "permutation", 
                               tuneGrid = data.frame(mtry = 50, min.node.size = 5, splitrule = "variance"))

# Select the top 100 important features
FeatureImportance <- sort(Differentiation.model$finalModel$variable.importance, decreasing = T)
SelectedFeatures <- names(FeatureImportance[1:100])
##   [1] "Aldoc"    "Stmn1"    "Slc1a2"   "Mt1"      "Sox11"    "Rtn1"     "Neurod1"  "Sox4"     "Tmem47"  
##  [10] "Basp1"    "Tmsb10"   "Stmn2"    "Tmsb4x"   "Nfib"     "Tubb5"    "Marcksl1" "Tubb3"    "Sparc"   
##  [19] "Calb2"    "Nfix"     "Marcks"   "Mt2"      "Dpysl3"   "Thrsp"    "Fxyd6"    "Myt1l"    "Ptma"    
##  [28] "Fabp7"    "Glul"     "Hnrnpa1"  "Zfp36l1"  "Kif5c"    "Ank3"     "Scrg1"    "Ass1"     "Mllt11"  
##  [37] "Gja1"     "Gfap"     "AW047730" "Ndufc2"   "Auts2"    "Cpe"      "Rnf165"   "Elavl3"   "Car2"    
##  [46] "Ly6h"     "Hmgn2"    "Kcnj10"   "Rplp0"    "Cxcl14"   "S100a6"   "Ppp1r14b" "Rpsa"     "Gas1"    
##  [55] "Kcnk1"    "Sept3"    "Hmgb3"    "Ier2"     "Elavl2"   "Baalc"    "Luzp2"    "Timp3"    "Rgcc"    
##  [64] "Gm10561"  "Hmgb2"    "Egr1"     "Eomes"    "Rps19"    "Draxin"   "Sorbs2"   "Sox1ot"   "H2afz"   
##  [73] "Nrep"     "Rpl29"    "Vim"      "Rpl32"    "Chchd10"  "Tril"     "Hes5"     "Itih3"    "Rab3a"   
##  [82] "Ccnd2"    "Aqp4"     "Mex3b"    "Gpc5"     "Zbtb18"   "Enc1"     "Stmn3"    "Csrp1"    "Gm6145"  
##  [91] "Rps2"     "Serpine2" "Stmn4"    "Lix1"     "Notch2"   "Eef1b2"   "Prox1"    "Ier5"     "Lmnb1"   
## [100] "Rps18"

# Re-train the model with the restricted features
Differentiation.model <- train(PseudoDifferentiation~., 
                               data = model.data[,colnames(model.data) %in% c(SelectedFeatures, "PseudoDifferentiation")], 
                               method = "ranger",
                               trControl = trctrl,
                               importance = "permutation", 
                               tuneGrid = data.frame(mtry = 50, min.node.size = 5, splitrule = "variance"))

# Plot Observed vs OOB predicted values from the model
Correlation <- cbind(model.data$PseudoDifferentiation, Differentiation.model$finalModel$predictions)
Correlation <- as.data.frame(Correlation)
ggplot(Correlation, aes(model.data$PseudoDifferentiation, Differentiation.model$finalModel$predictions)) + 
  geom_point(color = "#993404", alpha = 0.9, size = 0.8) + 
  geom_smooth(method = "lm", color = "black") + 
  ggtitle("") + xlab("Observed Pseudo-differentiation") + ylab("Predicted Pseudo-differentiation") + 
  theme_bw()
ggsave("Correlation_Observed_Predicted_Model.png", width = 4, height = 4)
ggsave("Correlation_Observed_Predicted_Model.pdf", width = 4, height = 4)

Model1 <- lm(Differentiation.model$finalModel$predictions ~ model.data$PseudoDifferentiation, data = Correlation)
summary(Model1)
round(summary(Model1)$adj.r.squared, digits = 2)
## [1] 0.99

######################################################################################################
# Predict differentiation score of qNSC
######################################################################################################
# Predict differentiation score of qNSC
qNSC.df <- as.data.frame(t(seurat_qNSC@assays[["RNA"]]@scale.data[rownames(seurat_qNSC@assays[["RNA"]]@scale.data) %in% 
                                                                    SelectedFeatures, ]))
qNSC.pred.maturation <- predict(Differentiation.model, qNSC.df)
seurat_qNSC$PseudoDifferentiation.score <- qNSC.pred.maturation

qNSC.plot.df <- seurat_qNSC@meta.data %>% 
  dplyr::select(Age, PseudoDifferentiation.score)
qNSC.plot.df$Age <- factor(x = qNSC.plot.df$Age, levels = c("Young", "Middle", "Old"))

# Plot differentiation score of qNSC by Age
ggplot(qNSC.plot.df, aes(x = Age, y = PseudoDifferentiation.score, fill = Age)) +
  geom_boxplot() +
  ## geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = c("#fcc5c0", "#f768a1", "#ae017e")) + 
  ggtitle("qNSC") + 
  xlab("") + 
  ylab("Differentiation score") + 
  theme_classic() + 
  NoLegend()
ggsave("Pseudo_differentiation_Prediction_qNSC.png", width = 4, height = 4)
ggsave("Pseudo_differentiation_Prediction_qNSC.pdf", width = 4, height = 4)









##== aNSPC ==##
##== With cell-cycle gene regress out ==##
##== Data process ==##
######################################################################################################
# Regress out cell-cycle genes (97 genes)
######################################################################################################
# Cell-cycle genes
s_genes <- c("Mcm5", "Pcna", "Tym5", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", 
             "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rap2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", 
             "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", 
             "Dscc1", "Blm", " Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
g2m_genes <- c("Hmgb2", "Ddk1","Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", 
               "Mki67", "Tmpo", " Cenpk", "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", 
               "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "kif20b", "Hjurp", "Cdca3", "Hn1", "Cdc20", "Ttk", 
               "Cdc25c", "kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", 
               "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")

# Regress out cell-cycle genes
exprMat <- as.matrix(seurat_postmitotic@assays$RNA@counts)
exprMat <- exprMat[!rownames(exprMat) %in% c(s_genes, g2m_genes), ]
dim(exprMat)
meta_postmitotic <- as.data.frame(seurat_postmitotic@meta.data)

# Re-clustering
seurat_postmitotic <- CreateSeuratObject(counts = exprMat)
DefaultAssay(seurat_postmitotic) <- "RNA"
seurat_postmitotic <- NormalizeData(seurat_postmitotic, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
seurat_postmitotic <- RunUMAP(seurat_postmitotic, reduction = "pca", dims = 1:15) %>% 
  RunTSNE(reduction = "pca", dims = 1:15) %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = 1.0)
seurat_postmitotic <- AddMetaData(seurat_postmitotic, metadata = meta_postmitotic)
Idents(seurat_postmitotic) <- "Celltype"
seurat_postmitotic$Celltype <- seurat_postmitotic@active.ident
DimPlot(seurat_postmitotic, reduction = "umap", cols = c("#fe9929", "#3690c0"), pt.size = 2)
ggsave("UMAP_Neurogenic_Postmitotic.png", width = 5.5, height = 4)
ggsave("UMAP_Neurogenic_Postmitotic.pdf", width = 5.5, height = 4)

FeaturePlot(seurat_postmitotic, features = c("Ascl1", "Sox2", "Dcx", "Neurod1"), 
            cols = c(rgb(0.7, 0.7, 0.7, alpha = 0.2), "darkred"), 
            ncol = 2) 
ggsave("Featureplot_Neurogenic_Postmitotic.png", width = 9, height = 8)
ggsave("Featureplot_Neurogenic_Postmitotic.pdf", width = 9, height = 8)

######################################################################################################
# Getting Started with Monocle
######################################################################################################
# Importing data from Seurat object
data <- as(as.matrix(seurat_postmitotic@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = seurat_postmitotic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

# Create CellDataSet object
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, 
                             relative_expr = T, 
                             cores = detectCores())

# Filtering low-quality cells
mycds <- detectGenes(mycds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(mycds), 
                                    num_cells_expressed >= 10)) 

######################################################################################################
# Constructing single cell trajectories
######################################################################################################
# Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(mycds[expressed_genes, ], 
                                      fullModelFormulaStr = "~Celltype", 
                                      cores = detectCores()) 
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) # 166 genes

mycds <- setOrderingFilter(mycds, ordering_genes)  
plot_ordering_genes(mycds)

# Trajectory step 2: reduce data dimensionality
mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree")

# Trajectory step 3: order cells along the trajectory
mycds <- orderCells(mycds, reverse = T)

######################################################################################################
# Differential Expression Analysis
######################################################################################################
# Clustering Genes by Pseudotemporal Expression Pattern
diff_test_res <- differentialGeneTest(mycds[ordering_genes, ],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01)) # 116 genes
p <- plot_pseudotime_heatmap(mycds[sig_gene_names, ],
                             num_clusters = 4,
                             cores = detectCores(),
                             show_rownames = T, 
                             hmcols = viridis(1000, option = "inferno", alpha = 1, begin = 0), 
                             return_heatmap = T)
ggsave("Monocle_Heatmap_Postmitotic.png", p, width = 8, height = 10)
ggsave("Monocle_Heatmap_Postmitotic.pdf", p, width = 8, height = 10)

# Extract pseudotemporally ordered genes in qNSCs
t <- as.data.frame(cutree(p$tree_row, k = 4))
colnames(t) <- "Cluster"
t$Cluster <- recode(t$Cluster, 
                    "1" = "module_1", 
                    "2" = "module_2", 
                    "3" = "module_3", 
                    "4" = "module_4")
t$Gene <- rownames(t)
module_1 <- t %>% 
  filter(Cluster == "module_1")
module_1 <- rownames(module_1)
module_2 <- t %>% 
  filter(Cluster == "module_2")
module_2 <- rownames(module_2)
module_3 <- t %>% 
  filter(Cluster == "module_3")
module_3 <- rownames(module_3)
module_4 <- t %>% 
  filter(Cluster == "module_4")
module_4 <- rownames(module_4)
write.csv(t, "pseudotime.Postmitotic.All.csv") # Used for GSEA

# Asign module score
DefaultAssay(seurat_postmitotic) <- "RNA"
seurat_postmitotic <- AddModuleScore(seurat_postmitotic, 
                                     features = list(module_1), 
                                     name = "module_1", 
                                     assay = "RNA")
seurat_postmitotic <- AddModuleScore(seurat_postmitotic, 
                                     features = list(module_2), 
                                     name = "module_2", 
                                     assay = "RNA")
seurat_postmitotic <- AddModuleScore(seurat_postmitotic, 
                                     features = list(module_3), 
                                     name = "module_3", 
                                     assay = "RNA")
seurat_postmitotic <- AddModuleScore(seurat_postmitotic, 
                                     features = list(module_4), 
                                     name = "module_4", 
                                     assay = "RNA")
seurat_postmitotic$module_1 <- seurat_postmitotic$module_11
seurat_postmitotic$module_2 <- seurat_postmitotic$module_21
seurat_postmitotic$module_3 <- seurat_postmitotic$module_31
seurat_postmitotic$module_4 <- seurat_postmitotic$module_41
seurat_postmitotic$module_11 <- NULL
seurat_postmitotic$module_21 <- NULL
seurat_postmitotic$module_31 <- NULL
seurat_postmitotic$module_41 <- NULL
Idents(seurat_postmitotic) <- "Celltype"
VlnPlot(seurat_postmitotic, features = c("module_1", "module_2", "module_3", "module_4"), 
        cols = c("#993404", "#fe9929", "#3690c0"), pt.size = 0, ncol = 4)
FeaturePlot(seurat_postmitotic, features = c("module_1", "module_2", "module_3", "module_4"), 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90", pt.size = 1.2, 
            cols = c(rgb(0.7, 0.7, 0.7, alpha = 0.2), "darkred"), ncol = 4)

