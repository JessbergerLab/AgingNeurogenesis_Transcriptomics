# Load data
Hippo_All <- readRDS("Data/SeuratObject/Visium/Processed/Hippo_All.rds")

######################################################################################################
# Module score for IFNg
######################################################################################################
# Create module score
library(msigdbr)
gene_set <- msigdbr(species = "Mus musculus", category = "H") %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(gene_symbol) %>%
  unique()
Hallmark_IFNg <- gene_set
Hippo_All <- AddModuleScore(Hippo_All, 
                            features = list(Hallmark_IFNg), 
                            name = "Hallmark_IFNg", 
                            assay = "SCT")
Hippo_All$Hallmark_IFNg <- Hippo_All$Hallmark_IFNg1
Hippo_All$Hallmark_IFNg1 <- NULL

# Plotting module score
ggplot(Hippo_All@meta.data, aes(x = Hallmark_IFNg)) + 
  geom_histogram(bins = 50) + 
  geom_vline(aes(xintercept = 0), color  = "red", linetype = 2) + 
  theme_bw()
ggsave("Histogram_Visium_IFNg.png", width = 4.5, height = 4.5)
ggsave("Histogram_Visium_IFNg.pdf", width = 4.5, height = 4.5)

# Subset population with module score > 0
Hippo_All_IFNg_Hi <- subset(Hippo_All, subset = Hallmark_IFNg > 0)
Hippo_All_IFNg_Hi <- WhichCells(Hippo_All_IFNg_Hi)
Idents(Hippo_All) <- "Region"
Idents(Hippo_All, cells = Hippo_All_IFNg_Hi) <- "Edge"
Hippo_All$FuncRegion <- Hippo_All@active.ident
Idents(Hippo_All) <- "FuncRegion"
Hippo_IFNg_Hi <- subset(Hippo_All, idents = "Edge")
metadata_Hippo_IFNg_Hi <- data.frame(Hippo_IFNg_Hi@meta.data)

table(Hippo_All$FuncRegion, Hippo_All$Age)
##      Young Middle  Old
## Edge    61     60  299
table(Hippo_All$Age)
## Young Middle    Old 
##  2483   2402   5036

# Export metadata for STutility
metadata_Hippo <- data.frame(Hippo_All@meta.data)
metadata_Hippo$ID <- rownames(metadata_Hippo)
metadata_Hippo$ID <- gsub("10_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("11_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("12_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("13_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("14_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("15_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("16_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("1_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("2_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("3_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("4_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("5_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("6_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("7_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("8_", "", metadata_Hippo$ID)
metadata_Hippo$ID <- gsub("9_", "", metadata_Hippo$ID)
rownames(metadata_Hippo) <- make.unique(metadata_Hippo$ID)





##== STutility ==##
######################################################################################################
# O_1_A1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_1/A1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_1/A1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_1/A1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_1/A1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_1_A1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "9")
rownames(metadata_O_1_A1) <- make.unique(metadata_O_1_A1$ID)
names <- make.unique(metadata_O_1_A1$ID)
rownames(metadata_O_1_A1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_1_A1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_1_A1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_1_A1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_1_A1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_1_A1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_1_A1 <- subset(RS_O_1_A1, !(RS_O_1_A1 %in% Edge_O_1_A1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_A1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_1_A1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_1_A1_2 <- subset(RS_O_1_A1_2, !(RS_O_1_A1_2 %in% Edge_O_1_A1))
RS_O_1_A1_2 <- subset(RS_O_1_A1_2, !(RS_O_1_A1_2 %in% nbs_Edge_O_1_A1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_A1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
## RS2      RS1 nbs_Edge     Edge      CA1      CA3 
## 129      140      119       45      205       59

# Export metadata
metadata_O_1_A1_Edge <- data.frame(se@meta.data)
metadata_O_1_A1_Edge <- metadata_O_1_A1_Edge[, -(25:26)]
metadata_O_1_A1_Edge <- metadata_O_1_A1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_1_A1_Edge$Barcode <- paste0(9, "_", metadata_O_1_A1_Edge$ID)

######################################################################################################
# O_1_B1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_1/B1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_1/B1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_1/B1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_1/B1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_1_B1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "10")
names <- make.unique(metadata_O_1_B1$ID)
rownames(metadata_O_1_B1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_1_B1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_1_B1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_1_B1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_1_B1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_1_B1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_1_B1 <- subset(RS_O_1_B1, !(RS_O_1_B1 %in% Edge_O_1_B1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_B1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_1_B1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_1_B1_2 <- subset(RS_O_1_B1_2, !(RS_O_1_B1_2 %in% Edge_O_1_B1))
RS_O_1_B1_2 <- subset(RS_O_1_B1_2, !(RS_O_1_B1_2 %in% nbs_Edge_O_1_B1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_B1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3       ML      GCL    Hilus 
##  113      169      166       51       14       12       29        6        5        3 

# Export metadata
metadata_O_1_B1_Edge <- data.frame(se@meta.data)
metadata_O_1_B1_Edge <- metadata_O_1_B1_Edge[, -(25:26)]
metadata_O_1_B1_Edge <- metadata_O_1_B1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_1_B1_Edge$Barcode <- paste0(10, "_", metadata_O_1_B1_Edge$ID)

######################################################################################################
# O_1_C1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_1/C1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_1/C1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_1/C1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_1/C1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_1_C1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "11")
names <- make.unique(metadata_O_1_C1$ID)
rownames(metadata_O_1_C1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_1_C1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_1_C1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_1_C1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_1_C1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_1_C1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_1_C1 <- subset(RS_O_1_C1, !(RS_O_1_C1 %in% Edge_O_1_C1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_C1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_1_C1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_1_C1_2 <- subset(RS_O_1_C1_2, !(RS_O_1_C1_2 %in% Edge_O_1_C1))
RS_O_1_C1_2 <- subset(RS_O_1_C1_2, !(RS_O_1_C1_2 %in% nbs_Edge_O_1_C1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_C1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3 
##  114      125      101       31       47       11      118

# Export metadata
metadata_O_1_C1_Edge <- data.frame(se@meta.data)
metadata_O_1_C1_Edge <- metadata_O_1_C1_Edge[, -(25:26)]
metadata_O_1_C1_Edge <- metadata_O_1_C1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_1_C1_Edge$Barcode <- paste0(11, "_", metadata_O_1_C1_Edge$ID)

######################################################################################################
# O_1_D1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_1/D1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_1/D1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_1/D1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_1/D1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_1_D1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "12")
names <- make.unique(metadata_O_1_D1$ID)
rownames(metadata_O_1_D1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_1_D1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_1_D1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_1_D1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_1_D1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_1_D1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_1_D1 <- subset(RS_O_1_D1, !(RS_O_1_D1 %in% Edge_O_1_D1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_D1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_1_D1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_1_D1_2 <- subset(RS_O_1_D1_2, !(RS_O_1_D1_2 %in% Edge_O_1_D1))
RS_O_1_D1_2 <- subset(RS_O_1_D1_2, !(RS_O_1_D1_2 %in% nbs_Edge_O_1_D1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_1_D1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3       ML      GCL    Hilus 
##  111      108       81       26      120       16       80       41       21        4 

# Export metadata
metadata_O_1_D1_Edge <- data.frame(se@meta.data)
metadata_O_1_D1_Edge <- metadata_O_1_D1_Edge[, -(25:26)]
metadata_O_1_D1_Edge <- metadata_O_1_D1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_1_D1_Edge$Barcode <- paste0(12, "_", metadata_O_1_D1_Edge$ID)

######################################################################################################
# O_2_A1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_2/A1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_2/A1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_2/A1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_2/A1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_2_A1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "13")
rownames(metadata_O_2_A1) <- make.unique(metadata_O_2_A1$ID)
names <- make.unique(metadata_O_2_A1$ID)
rownames(metadata_O_2_A1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_2_A1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_2_A1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_2_A1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_2_A1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_2_A1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_2_A1 <- subset(RS_O_2_A1, !(RS_O_2_A1 %in% Edge_O_2_A1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_A1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_2_A1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_2_A1_2 <- subset(RS_O_2_A1_2, !(RS_O_2_A1_2 %in% Edge_O_2_A1))
RS_O_2_A1_2 <- subset(RS_O_2_A1_2, !(RS_O_2_A1_2 %in% nbs_Edge_O_2_A1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_A1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA3       ML      GCL    Hilus 
##  142      161      128       44      139       37       12        9        9

# Export metadata
metadata_O_2_A1_Edge <- data.frame(se@meta.data)
metadata_O_2_A1_Edge <- metadata_O_2_A1_Edge[, -(25:26)]
metadata_O_2_A1_Edge <- metadata_O_2_A1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_2_A1_Edge$Barcode <- paste0(13, "_", metadata_O_2_A1_Edge$ID)

######################################################################################################
# O_2_B1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_2/B1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_2/B1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_2/B1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_2/B1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_2_B1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "14")
names <- make.unique(metadata_O_2_B1$ID)
rownames(metadata_O_2_B1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_2_B1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_2_B1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_2_B1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_2_B1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_2_B1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_2_B1 <- subset(RS_O_2_B1, !(RS_O_2_B1 %in% Edge_O_2_B1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_B1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_2_B1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_2_B1_2 <- subset(RS_O_2_B1_2, !(RS_O_2_B1_2 %in% Edge_O_2_B1))
RS_O_2_B1_2 <- subset(RS_O_2_B1_2, !(RS_O_2_B1_2 %in% nbs_Edge_O_2_B1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_B1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3       ML      GCL    Hilus 
##  137      198      199       56       40       12       38        5        3        3  

# Export metadata
metadata_O_2_B1_Edge <- data.frame(se@meta.data)
metadata_O_2_B1_Edge <- metadata_O_2_B1_Edge[, -(25:26)]
metadata_O_2_B1_Edge <- metadata_O_2_B1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_2_B1_Edge$Barcode <- paste0(14, "_", metadata_O_2_B1_Edge$ID)

######################################################################################################
# O_2_C1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_2/C1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_2/C1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_2/C1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_2/C1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_2_C1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "15")
names <- make.unique(metadata_O_2_C1$ID)
rownames(metadata_O_2_C1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_2_C1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_2_C1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_2_C1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_2_C1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_2_C1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_2_C1 <- subset(RS_O_2_C1, !(RS_O_2_C1 %in% Edge_O_2_C1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_C1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_2_C1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_2_C1_2 <- subset(RS_O_2_C1_2, !(RS_O_2_C1_2 %in% Edge_O_2_C1))
RS_O_2_C1_2 <- subset(RS_O_2_C1_2, !(RS_O_2_C1_2 %in% nbs_Edge_O_2_C1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_C1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
##  RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3       ML      GCL 
##  129      118       74       20      164       19       38        6        1

# Export metadata
metadata_O_2_C1_Edge <- data.frame(se@meta.data)
metadata_O_2_C1_Edge <- metadata_O_2_C1_Edge[, -(25:26)]
metadata_O_2_C1_Edge <- metadata_O_2_C1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_2_C1_Edge$Barcode <- paste0(15, "_", metadata_O_2_C1_Edge$ID)

######################################################################################################
# O_2_D1
######################################################################################################
# Load data
df <- data.frame(samples  =c("Data/Raw/Visium/O_2/D1/filtered_feature_bc_matrix.h5"),
                 spotfiles = c("Data/Raw/Visium/O_2/D1/spatial/tissue_positions_list.csv"),
                 imgs  = c("Data/Raw/Visium/O_2/D1/spatial/tissue_hires_image.png"),
                 json = c("Data/Raw/Visium/O_2/D1/spatial/scalefactors_json.json"))
seurat <- InputFromTable(infotable = df, 
                         min.gene.count = 0, 
                         min.gene.spots = 0,
                         min.spot.count = 0,
                         platform =  "Visium")
## Load metadata
metadata_O_2_D1 <- metadata_Hippo %>% 
  dplyr::filter(GEMgroup == "16")
names <- make.unique(metadata_O_2_D1$ID)
rownames(metadata_O_2_D1) <- paste0(names, "_", 1)
seurat <- AddMetaData(seurat, metadata = metadata_O_2_D1)
## Subset
head(seurat@meta.data)
Idents(seurat) <- "FuncRegion"
se <- subset(seurat, idents = c("Edge", "ML", "GCL", "Hilus", "CA1", "CA2", "CA3"))

# Load image
se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Neighborhood analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "FuncRegion")
Edge_O_2_D1 <- WhichCells(se, idents = "Edge")
se <- SetIdent(se, value = "nbs_Edge")
nbs_Edge_O_2_D1 <- WhichCells(se, idents = "nbs_Edge")
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = nbs_Edge_O_2_D1) <- "nbs_Edge" 
se$FuncRegion <- se@active.ident

# Remote spot1 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "nbs_Edge", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_nbs_Edge")
RS_O_2_D1 <- WhichCells(se, idents = "nbs_nbs_Edge")
RS_O_2_D1 <- subset(RS_O_2_D1, !(RS_O_2_D1 %in% Edge_O_2_D1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_D1) <- "RS1" 
se$FuncRegion <- se@active.ident

# Remote spot2 analysis
se <- SetIdent(se, value = "FuncRegion")
se <- RegionNeighbours(se, id = "RS1", verbose = TRUE)
## Extract spots
se <- SetIdent(se, value = "nbs_RS1")
RS_O_2_D1_2 <- WhichCells(se, idents = "nbs_RS1")
RS_O_2_D1_2 <- subset(RS_O_2_D1_2, !(RS_O_2_D1_2 %in% Edge_O_2_D1))
RS_O_2_D1_2 <- subset(RS_O_2_D1_2, !(RS_O_2_D1_2 %in% nbs_Edge_O_2_D1))
se <- SetIdent(se, value = "FuncRegion")
Idents(se, cells = RS_O_2_D1_2) <- "RS2" 
se$FuncRegion <- se@active.ident
table(se$FuncRegion)
## RS2      RS1 nbs_Edge     Edge      CA1      CA2      CA3       ML      GCL    Hilus 
## 130      144       99       26      127       71       42       20       11        5 

# Export metadata
metadata_O_2_D1_Edge <- data.frame(se@meta.data)
metadata_O_2_D1_Edge <- metadata_O_2_D1_Edge[, -(25:26)]
metadata_O_2_D1_Edge <- metadata_O_2_D1_Edge %>% 
  filter(FuncRegion != "CA1") %>% 
  filter(FuncRegion != "CA2") %>% 
  filter(FuncRegion != "CA3") %>% 
  filter(FuncRegion != "ML") %>% 
  filter(FuncRegion != "GCL") %>% 
  filter(FuncRegion != "Hilus")
metadata_O_2_D1_Edge$Barcode <- paste0(16, "_", metadata_O_2_D1_Edge$ID)

######################################################################################################
# Merge all metadata with Edge-Neighbor information
######################################################################################################
metadata_Edge_Old <- rbind(metadata_O_1_A1_Edge, metadata_O_1_B1_Edge, metadata_O_1_C1_Edge, metadata_O_1_D1_Edge, 
                           metadata_O_2_A1_Edge, metadata_O_2_B1_Edge, metadata_O_2_C1_Edge, metadata_O_2_D1_Edge)
metadata_Edge_Old <- metadata_Edge_Old[, c(10, 23, 24, 26)]
rownames(metadata_Edge_Old) <- make.unique(metadata_Edge_Old$Barcode)

# Save the work space
write.csv(metadata_Edge_Old, "metadata_Edge_Old_IS_NNS_RS_2.csv")






##== Inflammatory gradient ==##
# Load data
Hippo_All <- readRDS("Data/SeuratObject/Visium/Processed/Hippo_All.rds")
metadata_Edge <- read.csv("./Data/DEG/Visium/metadata_Edge_Old_IS_NNS_RS_2.csv")
rownames(metadata_Edge) <- metadata_Edge[, 1]
metadata_Edge <- metadata_Edge[, -1]

seurat <- AddMetaData(Hippo_All, metadata = metadata_Edge)
Idents(seurat) <- "FuncRegion"
seurat <- subset(seurat, idents = c("Edge", "nbs_Edge", "RS1", "RS2"))
seurat$NBS <- recode(seurat$FuncRegion, 
                     "Edge" = "IS", 
                     "nbs_Edge" = "NNS", 
                     "RS1" = "RS1", 
                     "RS2" = "RS2")
seurat$Group <- paste0(seurat$GEMgroup, "_", seurat$NBS)
seurat$Group <- recode(seurat$Group, 
                       "9_IS" = "1", 
                       "9_NNS" = "2", 
                       "9_RS1" = "3", 
                       "9_RS2" = "4", 
                       "10_IS" = "5", 
                       "10_NNS" = "6", 
                       "10_RS1" = "7", 
                       "10_RS2" = "8", 
                       "11_IS" = "9", 
                       "11_NNS" = "10", 
                       "11_RS1" = "11", 
                       "11_RS2" = "12", 
                       "12_IS" = "13", 
                       "12_NNS" = "14", 
                       "12_RS1" = "15", 
                       "12_RS2" = "16", 
                       "13_IS" = "17", 
                       "13_NNS" = "18", 
                       "13_RS1" = "19", 
                       "13_RS2" = "20", 
                       "14_IS" = "21", 
                       "14_NNS" = "22", 
                       "14_RS1" = "23", 
                       "14_RS2" = "24", 
                       "15_IS" = "25", 
                       "15_NNS" = "26", 
                       "15_RS1" = "27", 
                       "15_RS2" = "28", 
                       "16_IS" = "29", 
                       "16_NNS" = "30", 
                       "16_RS1" = "31", 
                       "16_RS2" = "32")
seurat$FuncRegion <- recode(seurat$FuncRegion, 
                            "Edge" = "IS", 
                            "nbs_Edge" = "NNS", 
                            "RS1" = "ENS1", 
                            "RS2" = "ENS2")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$Spatial@counts
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat$Age)
metadata$group_id <- factor(seurat$FuncRegion)
metadata$sample_id <- factor(seurat$Group)

######################################################################################################
# Create single cell experiment object
######################################################################################################
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

######################################################################################################
# Acquiring necessary metrics for aggregation across cells in a sample
######################################################################################################
# Named vector of cluster names (Region)
cids <- purrr::set_names(levels(sce$cluster_id))
## Total number of clusters
nk <- length(cids)
nk

# Named vector of sample names (GEMgroup)
sids <- purrr::set_names(levels(sce$sample_id))
## Total number of samples 
ns <- length(sids)
ns

# Generate sample level metadata
## Determine the number of cells per sample
table(sce$sample_id)
## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei

######################################################################################################
# Count aggregation to sample level
######################################################################################################
# Aggregate the counts per sample_id and cluster_id
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
class(pb)
dim(pb)
pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)

######################################################################################################
# Save the work space
######################################################################################################
save.image(file = 'DESeq2_Visium_IFNg_RemoteSpot2.RData')









load('Data/SeuratObject/Environment/DESeq2/DESeq2_Visium_IFNg_RemoteSpot2.RData')
##==  Differential gene expression with DESeq2 ==##
######################################################################################################
# Sample-level metadata
######################################################################################################
# Get sample names for each of the cell type clusters
# Prepare data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}
de_samples <- map(1:length(cids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples
samples_list <- map(1:length(cids), get_sample_ids)
get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(cids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 

metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 
metadata <- as.data.frame(metadata)

######################################################################################################
# Subsetting Old
######################################################################################################
# Generate vector of cluster IDs
clusters <- unique(metadata$cluster_id)
clusters
clusters[1]

# Subset the metadata to only the DG
cluster_metadata <- metadata[which(metadata$cluster_id == "Old"), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the DG
counts <- pb[["Old"]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))

######################################################################################################
# Create DESeq2 object
######################################################################################################
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "group_id") + 
  scale_color_manual(breaks = c("IS", "NNS", "ENS1", "ENS2"), 
                     values = c("#e31a1c", "#fd8d3c", "#fed976", "#ffffcc")) + 
  theme_bw() + 
  coord_fixed(ratio = 3)
DESeq2::plotPCA(rld, intgroup = "group_id") + 
  scale_color_manual(breaks = c("IS", "NNS", "ENS1", "ENS2"), 
                     values = c("#e31a1c", "#fd8d3c", "#fed976", "#fed976")) + 
  theme_bw() + 
  coord_fixed(ratio = 3)
ggsave("DESeq2_PCA_Visium_InflammatoryGradient2.png", width = 5, height = 5)
ggsave("DESeq2_PCA_Visium_InflammatoryGradient2.pdf", width = 5, height = 5)

######################################################################################################
# Hierarchical clustering
######################################################################################################
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
ann_colors = c("#e31a1c", "#fd8d3c", "#fed976", "#ffeda0")

pheatmap(rld_cor, 
         annotation = cluster_metadata[, c("group_id"), drop = F], 
         show_colnames = F, show_rawnames = F, cellwidth = 18, cellheight = 18, 
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100), 
         filename = "Correlation_Visium_InflammatoryGradient2.pdf",
         width = 16, height = 16)

######################################################################################################
# Run DESeq2 differential expression analysis
######################################################################################################
dds <- DESeq(dds)
# Plot dispersion estimates
pdf("DispersionPlot_Visium_InflammatoryGradient2.pdf", width = 5, height = 4)
plotDispEsts(dds)
dev.off()




##== Visualization ==##
######################################################################################################
# VolcanoPlot 
######################################################################################################
# InflammatoryEdge
df <- results(dds, contrast = c('group_id', 'IS', 'NNS'))
head(df, 3)
df_selected <- as.data.frame(df) %>% 
  drop_na()
df_selected$gene <- rownames(df_selected)
df_selected$Change <- ifelse(df_selected$padj < 0.05,
                             ifelse(abs(df_selected$log2FoldChange) > log2(1.1),
                                    ifelse(df_selected$log2FoldChange < -log2(1.1), 'Down', 'Up'), 'Stable'), 'Stable')
## Check
table(df_selected$Change)
## Down Stable     Up 
## 1326  12138   2357 
## Add z-score based on two sided null hypothesis
df_selected$z <- p.to.Z(df_selected$pvalue) * sign(df_selected$log2FoldChange) 
df_selected$z.adj <- p.to.Z(df_selected$padj) * sign(df_selected$log2FoldChange)

top50 <- df_selected[order(abs(df_selected$z), decreasing = T), ]
top50 <- top50 %>% 
  filter(Change != "Stable") 
top50 <- top50[1:50, ]
## Volcano plot
ggplot(df_selected, aes(x = log2FoldChange, y = -log10(padj), fill = "Change")) +
  geom_point(aes(colour = Change)) + 
  scale_color_manual(values = c("#fd8d3c", "lightgrey", "#e31a1c")) + 
  ggtitle("Inflammatory Spots vs. Nearest Neighbor Spots") +
  xlab(log[2]~FC) + 
  ylab(-log[10]~p_adj_val) + 
  geom_hline(aes(yintercept = -log10(0.05)), color  = "darkgrey", linetype = "dashed") +
  geom_vline(aes(xintercept = log2(1.1)), color  = "darkgrey", linetype = "dashed") + 
  geom_vline(aes(xintercept = -log2(1.1)), color  = "darkgrey", linetype = "dashed") + 
  geom_text_repel(data = head(top50, 50), aes(label = gene)) + 
  theme_classic() + 
  NoLegend()
ggsave("DESeq2_Visium_Volcanoplot_InflammatoryEdge.png", width = 9, height = 8)
ggsave("DESeq2_Visium_Volcanoplot_InflammatoryEdge.pdf", width = 9, height = 8)

