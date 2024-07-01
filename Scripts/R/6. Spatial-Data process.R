##== Data process ==##
######################################################################################################
# Load library and data
######################################################################################################
# Load library
rm(list=ls())
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
# Load data
data_dir_Y_A1 <- "./Data/Raw/Visium/Y/A1/"
list.files(data_dir_Y_A1) # Should show filtered_feature_bc_matrix.h5
Y_A1 <- Load10X_Spatial(data.dir = data_dir_Y_A1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_Y_B1 <- "./Data/Raw/Visium/Y/B1/"
list.files(data_dir_Y_B1) # Should show filtered_feature_bc_matrix.h5
Y_B1 <- Load10X_Spatial(data.dir = data_dir_Y_B1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_Y_C1 <- "./Data/Raw/Visium/Y/C1/"
list.files(data_dir_Y_C1) # Should show filtered_feature_bc_matrix.h5
Y_C1 <- Load10X_Spatial(data.dir = data_dir_Y_C1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_Y_D1 <- "./Data/Raw/Visium/Y/D1/"
list.files(data_dir_Y_D1) # Should show filtered_feature_bc_matrix.h5
Y_D1 <- Load10X_Spatial(data.dir = data_dir_Y_D1, 
                        filename = "filtered_feature_bc_matrix.h5")

data_dir_M_A1 <- "./Data/Raw/Visium/M/A1/"
list.files(data_dir_M_A1) # Should show filtered_feature_bc_matrix.h5
M_A1 <- Load10X_Spatial(data.dir = data_dir_M_A1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_M_B1 <- "./Data/Raw/Visium/M/B1/"
list.files(data_dir_M_B1) # Should show filtered_feature_bc_matrix.h5
M_B1 <- Load10X_Spatial(data.dir = data_dir_M_B1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_M_C1 <- "./Data/Raw/Visium/M/C1/"
list.files(data_dir_M_C1) # Should show filtered_feature_bc_matrix.h5
M_C1 <- Load10X_Spatial(data.dir = data_dir_M_C1, 
                        filename = "filtered_feature_bc_matrix.h5")
data_dir_M_D1 <- "./Data/Raw/Visium/M/D1/"
list.files(data_dir_M_D1) # Should show filtered_feature_bc_matrix.h5
M_D1 <- Load10X_Spatial(data.dir = data_dir_M_D1, 
                        filename = "filtered_feature_bc_matrix.h5")

data_dir_O_1_A1 <- "./Data/Raw/Visium/O_1/A1/"
list.files(data_dir_O_1_A1) # Should show filtered_feature_bc_matrix.h5
O_1_A1 <- Load10X_Spatial(data.dir = data_dir_O_1_A1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_1_B1 <- "./Data/Raw/Visium/O_1/B1/"
list.files(data_dir_O_1_B1) # Should show filtered_feature_bc_matrix.h5
O_1_B1 <- Load10X_Spatial(data.dir = data_dir_O_1_B1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_1_C1 <- "./Data/Raw/Visium/O_1/C1/"
list.files(data_dir_O_1_C1) # Should show filtered_feature_bc_matrix.h5
O_1_C1 <- Load10X_Spatial(data.dir = data_dir_O_1_C1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_1_D1 <- "./Data/Raw/Visium/O_1/D1/"
list.files(data_dir_O_1_D1) # Should show filtered_feature_bc_matrix.h5
O_1_D1 <- Load10X_Spatial(data.dir = data_dir_O_1_D1, 
                          filename = "filtered_feature_bc_matrix.h5")

data_dir_O_2_A1 <- "./Data/Raw/Visium/O_2/A1/"
list.files(data_dir_O_2_A1) # Should show filtered_feature_bc_matrix.h5
O_2_A1 <- Load10X_Spatial(data.dir = data_dir_O_2_A1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_2_B1 <- "./Data/Raw/Visium/O_2/B1/"
list.files(data_dir_O_2_B1) # Should show filtered_feature_bc_matrix.h5
O_2_B1 <- Load10X_Spatial(data.dir = data_dir_O_2_B1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_2_C1 <- "./Data/Raw/Visium/O_2/C1/"
list.files(data_dir_O_2_C1) # Should show filtered_feature_bc_matrix.h5
O_2_C1 <- Load10X_Spatial(data.dir = data_dir_O_2_C1, 
                          filename = "filtered_feature_bc_matrix.h5")
data_dir_O_2_D1 <- "./Data/Raw/Visium/O_2/D1/"
list.files(data_dir_O_2_D1) # Should show filtered_feature_bc_matrix.h5
O_2_D1 <- Load10X_Spatial(data.dir = data_dir_O_2_D1, 
                          filename = "filtered_feature_bc_matrix.h5")

######################################################################################################
# Transfer the label of Artifacts
######################################################################################################
## Young
matadata_Y_A1 <- read.csv("./Data/Visium coordinates/Artifacts/Y_A1_Artifacts.csv")
names <- make.unique(matadata_Y_A1$Barcode)
rownames(matadata_Y_A1) <- names
colnames(matadata_Y_A1)[2] <- "Artifacts"
Y_A1 <- AddMetaData(Y_A1, metadata = matadata_Y_A1)
Y_A1$Barcode <- NULL
# Subset
Idents(Y_A1) <- "Artifacts"
Y_A1$Artifacts <- "NA"
Idents(Y_A1) <- "Artifacts"
Idents(Y_A1, cells = names) <- "Yes" 
Y_A1$Artifacts <- Y_A1@active.ident
table(Y_A1$Artifacts)
Y_A1_clean <- subset(Y_A1, idents = "Yes", invert = T)

matadata_Y_B1 <- read.csv("./Data/Visium coordinates/Artifacts/Y_B1_Artifacts.csv")
names <- make.unique(matadata_Y_B1$Barcode)
rownames(matadata_Y_B1) <- names
colnames(matadata_Y_B1)[2] <- "Artifacts"
Y_B1 <- AddMetaData(Y_B1, metadata = matadata_Y_B1)
Y_B1$Barcode <- NULL
# Subset
Idents(Y_B1) <- "Artifacts"
Y_B1$Artifacts <- "NA"
Idents(Y_B1) <- "Artifacts"
Idents(Y_B1, cells = names) <- "Yes" 
Y_B1$Artifacts <- Y_B1@active.ident
table(Y_B1$Artifacts)
Y_B1_clean <- subset(Y_B1, idents = "Yes", invert = T)

matadata_Y_C1 <- read.csv("./Data/Visium coordinates/Artifacts/Y_C1_Artifacts.csv")
names <- make.unique(matadata_Y_C1$Barcode)
rownames(matadata_Y_C1) <- names
colnames(matadata_Y_C1)[2] <- "Artifacts"
Y_C1 <- AddMetaData(Y_C1, metadata = matadata_Y_C1)
Y_C1$Barcode <- NULL
# Subset
Idents(Y_C1) <- "Artifacts"
Y_C1$Artifacts <- "NA"
Idents(Y_C1) <- "Artifacts"
Idents(Y_C1, cells = names) <- "Yes" 
Y_C1$Artifacts <- Y_C1@active.ident
table(Y_C1$Artifacts)
Y_C1_clean <- subset(Y_C1, idents = "Yes", invert = T)

matadata_Y_D1 <- read.csv("./Data/Visium coordinates/Artifacts/Y_D1_Artifacts.csv")
names <- make.unique(matadata_Y_D1$Barcode)
rownames(matadata_Y_D1) <- names
colnames(matadata_Y_D1)[2] <- "Artifacts"
Y_D1 <- AddMetaData(Y_D1, metadata = matadata_Y_D1)
Y_D1$Barcode <- NULL
# Subset
Idents(Y_D1) <- "Artifacts"
Y_D1$Artifacts <- "NA"
Idents(Y_D1) <- "Artifacts"
Idents(Y_D1, cells = names) <- "Yes" 
Y_D1$Artifacts <- Y_D1@active.ident
table(Y_D1$Artifacts)
Y_D1_clean <- subset(Y_D1, idents = "Yes", invert = T)

## Middle-age
matadata_M_A1 <- read.csv("./Data/Visium coordinates/Artifacts/M_A1_Artifacts.csv")
names <- make.unique(matadata_M_A1$Barcode)
rownames(matadata_M_A1) <- names
colnames(matadata_M_A1)[2] <- "Artifacts"
M_A1 <- AddMetaData(M_A1, metadata = matadata_M_A1)
M_A1$Barcode <- NULL
# Subset
Idents(M_A1) <- "Artifacts"
M_A1$Artifacts <- "NA"
Idents(M_A1) <- "Artifacts"
Idents(M_A1, cells = names) <- "Yes" 
M_A1$Artifacts <- M_A1@active.ident
table(M_A1$Artifacts)
M_A1_clean <- subset(M_A1, idents = "Yes", invert = T)

matadata_M_B1 <- read.csv("./Data/Visium coordinates/Artifacts/M_B1_Artifacts.csv")
names <- make.unique(matadata_M_B1$Barcode)
rownames(matadata_M_B1) <- names
colnames(matadata_M_B1)[2] <- "Artifacts"
M_B1 <- AddMetaData(M_B1, metadata = matadata_M_B1)
M_B1$Barcode <- NULL
# Subset
Idents(M_B1) <- "Artifacts"
M_B1$Artifacts <- "NA"
Idents(M_B1) <- "Artifacts"
Idents(M_B1, cells = names) <- "Yes" 
M_B1$Artifacts <- M_B1@active.ident
table(M_B1$Artifacts)
M_B1_clean <- subset(M_B1, idents = "Yes", invert = T)

matadata_M_C1 <- read.csv("./Data/Visium coordinates/Artifacts/M_C1_Artifacts.csv")
names <- make.unique(matadata_M_C1$Barcode)
rownames(matadata_M_C1) <- names
colnames(matadata_M_C1)[2] <- "Artifacts"
M_C1 <- AddMetaData(M_C1, metadata = matadata_M_C1)
M_C1$Barcode <- NULL
# Subset
Idents(M_C1) <- "Artifacts"
M_C1$Artifacts <- "NA"
Idents(M_C1) <- "Artifacts"
Idents(M_C1, cells = names) <- "Yes" 
M_C1$Artifacts <- M_C1@active.ident
table(M_C1$Artifacts)
M_C1_clean <- subset(M_C1, idents = "Yes", invert = T)

matadata_M_D1 <- read.csv("./Data/Visium coordinates/Artifacts/M_D1_Artifacts.csv")
names <- make.unique(matadata_M_D1$Barcode)
rownames(matadata_M_D1) <- names
colnames(matadata_M_D1)[2] <- "Artifacts"
M_D1 <- AddMetaData(M_D1, metadata = matadata_M_D1)
M_D1$Barcode <- NULL
# Subset
Idents(M_D1) <- "Artifacts"
M_D1$Artifacts <- "NA"
Idents(M_D1) <- "Artifacts"
Idents(M_D1, cells = names) <- "Yes" 
M_D1$Artifacts <- M_D1@active.ident
table(M_D1$Artifacts)
M_D1_clean <- subset(M_D1, idents = "Yes", invert = T)

## Old
matadata_O_1_A1 <- read.csv("./Data/Visium coordinates/Artifacts/O_1_A1_Artifacts.csv")
names <- make.unique(matadata_O_1_A1$Barcode)
rownames(matadata_O_1_A1) <- names
colnames(matadata_O_1_A1)[2] <- "Artifacts"
O_1_A1 <- AddMetaData(O_1_A1, metadata = matadata_O_1_A1)
O_1_A1$Barcode <- NULL
# Subset
Idents(O_1_A1) <- "Artifacts"
O_1_A1$Artifacts <- "NA"
Idents(O_1_A1) <- "Artifacts"
Idents(O_1_A1, cells = names) <- "Yes" 
O_1_A1$Artifacts <- O_1_A1@active.ident
table(O_1_A1$Artifacts)
O_1_A1_clean <- subset(O_1_A1, idents = "Yes", invert = T)

matadata_O_1_B1 <- read.csv("./Data/Visium coordinates/Artifacts/O_1_B1_Artifacts.csv")
names <- make.unique(matadata_O_1_B1$Barcode)
rownames(matadata_O_1_B1) <- names
colnames(matadata_O_1_B1)[2] <- "Artifacts"
O_1_B1 <- AddMetaData(O_1_B1, metadata = matadata_O_1_B1)
O_1_B1$Barcode <- NULL
# Subset
Idents(O_1_B1) <- "Artifacts"
O_1_B1$Artifacts <- "NA"
Idents(O_1_B1) <- "Artifacts"
Idents(O_1_B1, cells = names) <- "Yes" 
O_1_B1$Artifacts <- O_1_B1@active.ident
table(O_1_B1$Artifacts)
O_1_B1_clean <- subset(O_1_B1, idents = "Yes", invert = T)

matadata_O_1_C1 <- read.csv("./Data/Visium coordinates/Artifacts/O_1_C1_Artifacts.csv")
names <- make.unique(matadata_O_1_C1$Barcode)
rownames(matadata_O_1_C1) <- names
colnames(matadata_O_1_C1)[2] <- "Artifacts"
O_1_C1 <- AddMetaData(O_1_C1, metadata = matadata_O_1_C1)
O_1_C1$Barcode <- NULL
# Subset
Idents(O_1_C1) <- "Artifacts"
O_1_C1$Artifacts <- "NA"
Idents(O_1_C1) <- "Artifacts"
Idents(O_1_C1, cells = names) <- "Yes" 
O_1_C1$Artifacts <- O_1_C1@active.ident
table(O_1_C1$Artifacts)
O_1_C1_clean <- subset(O_1_C1, idents = "Yes", invert = T)

matadata_O_1_D1 <- read.csv("./Data/Visium coordinates/Artifacts/O_1_D1_Artifacts.csv")
names <- make.unique(matadata_O_1_D1$Barcode)
rownames(matadata_O_1_D1) <- names
colnames(matadata_O_1_D1)[2] <- "Artifacts"
O_1_D1 <- AddMetaData(O_1_D1, metadata = matadata_O_1_D1)
O_1_D1$Barcode <- NULL
# Subset
Idents(O_1_D1) <- "Artifacts"
O_1_D1$Artifacts <- "NA"
Idents(O_1_D1) <- "Artifacts"
Idents(O_1_D1, cells = names) <- "Yes" 
O_1_D1$Artifacts <- O_1_D1@active.ident
table(O_1_D1$Artifacts)
O_1_D1_clean <- subset(O_1_D1, idents = "Yes", invert = T)

matadata_O_2_A1 <- read.csv("./Data/Visium coordinates/Artifacts/O_2_A1_Artifacts.csv")
names <- make.unique(matadata_O_2_A1$Barcode)
rownames(matadata_O_2_A1) <- names
colnames(matadata_O_2_A1)[2] <- "Artifacts"
O_2_A1 <- AddMetaData(O_2_A1, metadata = matadata_O_2_A1)
O_2_A1$Barcode <- NULL
# Subset
Idents(O_2_A1) <- "Artifacts"
O_2_A1$Artifacts <- "NA"
Idents(O_2_A1) <- "Artifacts"
Idents(O_2_A1, cells = names) <- "Yes" 
O_2_A1$Artifacts <- O_2_A1@active.ident
table(O_2_A1$Artifacts)
O_2_A1_clean <- subset(O_2_A1, idents = "Yes", invert = T)

matadata_O_2_B1 <- read.csv("./Data/Visium coordinates/Artifacts/O_2_B1_Artifacts.csv")
names <- make.unique(matadata_O_2_B1$Barcode)
rownames(matadata_O_2_B1) <- names
colnames(matadata_O_2_B1)[2] <- "Artifacts"
O_2_B1 <- AddMetaData(O_2_B1, metadata = matadata_O_2_B1)
O_2_B1$Barcode <- NULL
# Subset
Idents(O_2_B1) <- "Artifacts"
O_2_B1$Artifacts <- "NA"
Idents(O_2_B1) <- "Artifacts"
Idents(O_2_B1, cells = names) <- "Yes" 
O_2_B1$Artifacts <- O_2_B1@active.ident
table(O_2_B1$Artifacts)
O_2_B1_clean <- subset(O_2_B1, idents = "Yes", invert = T)

matadata_O_2_C1 <- read.csv("./Data/Visium coordinates/Artifacts/O_2_C1_Artifacts.csv")
names <- make.unique(matadata_O_2_C1$Barcode)
rownames(matadata_O_2_C1) <- names
colnames(matadata_O_2_C1)[2] <- "Artifacts"
O_2_C1 <- AddMetaData(O_2_C1, metadata = matadata_O_2_C1)
O_2_C1$Barcode <- NULL
# Subset
Idents(O_2_C1) <- "Artifacts"
O_2_C1$Artifacts <- "NA"
Idents(O_2_C1) <- "Artifacts"
Idents(O_2_C1, cells = names) <- "Yes" 
O_2_C1$Artifacts <- O_2_C1@active.ident
table(O_2_C1$Artifacts)
O_2_C1_clean <- subset(O_2_C1, idents = "Yes", invert = T)

matadata_O_2_D1 <- read.csv("./Data/Visium coordinates/Artifacts/O_2_D1_Artifacts.csv")
names <- make.unique(matadata_O_2_D1$Barcode)
rownames(matadata_O_2_D1) <- names
colnames(matadata_O_2_D1)[2] <- "Artifacts"
O_2_D1 <- AddMetaData(O_2_D1, metadata = matadata_O_2_D1)
O_2_D1$Barcode <- NULL
# Subset
Idents(O_2_D1) <- "Artifacts"
O_2_D1$Artifacts <- "NA"
Idents(O_2_D1) <- "Artifacts"
Idents(O_2_D1, cells = names) <- "Yes" 
O_2_D1$Artifacts <- O_2_D1@active.ident
table(O_2_D1$Artifacts)
O_2_D1_clean <- subset(O_2_D1, idents = "Yes", invert = T)

######################################################################################################
# Perform SCT (recommended)
######################################################################################################
Y_A1_clean <- SCTransform(Y_A1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
Y_B1_clean <- SCTransform(Y_B1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
Y_C1_clean <- SCTransform(Y_C1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
Y_D1_clean <- SCTransform(Y_D1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
M_A1_clean <- SCTransform(M_A1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
M_B1_clean <- SCTransform(M_B1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
M_C1_clean <- SCTransform(M_C1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
M_D1_clean <- SCTransform(M_D1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_1_A1_clean <- SCTransform(O_1_A1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_1_B1_clean <- SCTransform(O_1_B1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_1_C1_clean <- SCTransform(O_1_C1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_1_D1_clean <- SCTransform(O_1_D1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_2_A1_clean <- SCTransform(O_2_A1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_2_B1_clean <- SCTransform(O_2_B1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_2_C1_clean <- SCTransform(O_2_C1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
O_2_D1_clean <- SCTransform(O_2_D1_clean, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
# Save results
saveRDS(Y_A1_clean, "Y_A1_clean.rds")
saveRDS(Y_B1_clean, "Y_B1_clean.rds")
saveRDS(Y_C1_clean, "Y_C1_clean.rds")
saveRDS(Y_D1_clean, "Y_D1_clean.rds")
saveRDS(M_A1_clean, "M_A1_clean.rds")
saveRDS(M_B1_clean, "M_B1_clean.rds")
saveRDS(M_C1_clean, "M_C1_clean.rds")
saveRDS(M_D1_clean, "M_D1_clean.rds")
saveRDS(O_1_A1_clean, "O_1_A1_clean.rds")
saveRDS(O_1_B1_clean, "O_1_B1_clean.rds")
saveRDS(O_1_C1_clean, "O_1_C1_clean.rds")
saveRDS(O_1_D1_clean, "O_1_D1_clean.rds")
saveRDS(O_2_A1_clean, "O_2_A1_clean.rds")
saveRDS(O_2_B1_clean, "O_2_B1_clean.rds")
saveRDS(O_2_C1_clean, "O_2_C1_clean.rds")
saveRDS(O_2_D1_clean, "O_2_D1_clean.rds")

######################################################################################################
# Integration of multislices
######################################################################################################
Whole_clean.merge <- merge(Y_A1_clean, y = c(Y_B1_clean, Y_C1_clean, Y_D1_clean, M_A1_clean, M_B1_clean, M_C1_clean, M_D1_clean, 
                                             O_1_A1_clean, O_1_B1_clean, O_1_C1_clean, O_1_D1_clean, O_2_A1_clean, O_2_B1_clean, O_2_C1_clean, O_2_D1_clean), 
                           add.cell.ids = c("1", "2", "3", "4", "5", "6", "7", "8", 
                                            "9", "10", "11", "12", "13", "14", "15", "16"), 
                           project = "WholeDG_clean")
Whole_clean.merge@meta.data$GEMgroup[grepl("1_", colnames(Whole_clean.merge))] <- "1"
Whole_clean.merge@meta.data$GEMgroup[grepl("2_", colnames(Whole_clean.merge))] <- "2"
Whole_clean.merge@meta.data$GEMgroup[grepl("3_", colnames(Whole_clean.merge))] <- "3"
Whole_clean.merge@meta.data$GEMgroup[grepl("4_", colnames(Whole_clean.merge))] <- "4"
Whole_clean.merge@meta.data$GEMgroup[grepl("5_", colnames(Whole_clean.merge))] <- "5"
Whole_clean.merge@meta.data$GEMgroup[grepl("6_", colnames(Whole_clean.merge))] <- "6"
Whole_clean.merge@meta.data$GEMgroup[grepl("7_", colnames(Whole_clean.merge))] <- "7"
Whole_clean.merge@meta.data$GEMgroup[grepl("8_", colnames(Whole_clean.merge))] <- "8"
Whole_clean.merge@meta.data$GEMgroup[grepl("9_", colnames(Whole_clean.merge))] <- "9"
Whole_clean.merge@meta.data$GEMgroup[grepl("10_", colnames(Whole_clean.merge))] <- "10"
Whole_clean.merge@meta.data$GEMgroup[grepl("11_", colnames(Whole_clean.merge))] <- "11"
Whole_clean.merge@meta.data$GEMgroup[grepl("12_", colnames(Whole_clean.merge))] <- "12"
Whole_clean.merge@meta.data$GEMgroup[grepl("13_", colnames(Whole_clean.merge))] <- "13"
Whole_clean.merge@meta.data$GEMgroup[grepl("14_", colnames(Whole_clean.merge))] <- "14"
Whole_clean.merge@meta.data$GEMgroup[grepl("15_", colnames(Whole_clean.merge))] <- "15"
Whole_clean.merge@meta.data$GEMgroup[grepl("16_", colnames(Whole_clean.merge))] <- "16"
my_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")
Whole_clean.merge@meta.data$GEMgroup <- factor(x = Whole_clean.merge@meta.data$GEMgroup, levels = my_levels)

# Add Position
Whole_clean.merge@meta.data$Position[grepl("1_|3_|5_|7_|9_|11_|13_|15_", colnames(Whole_clean.merge))] <- "Posterior"
Whole_clean.merge@meta.data$Position[grepl("2_|4_|6_|8_|10_|12_|14_|16_", colnames(Whole_clean.merge))] <- "Anterior"

# Add Age
Whole_clean.merge@meta.data$Age[grepl("1_|2_|3_|4_", colnames(Whole_clean.merge))] <- "Young"
Whole_clean.merge@meta.data$Age[grepl("5_|6_|7_|8_", colnames(Whole_clean.merge))] <- "Middle"
Whole_clean.merge@meta.data$Age[grepl("9_|10_|11_|12_|13_|14_|15_|16_", colnames(Whole_clean.merge))] <- "Old"

# Add MT & Ribo content 
Whole_clean.merge[["percentMito"]] <- PercentageFeatureSet(Whole_clean.merge, pattern = "^mt-")
Whole_clean.merge[["percentRibo"]] <- PercentageFeatureSet(Whole_clean.merge, pattern = "^Rps|^Rpl")

# Plotting basic information of Visium data
VlnPlot(Whole_clean.merge, features = "nCount_Spatial", group.by = "GEMgroup", pt.size = 0) + NoLegend()
ggsave("Vlnplot_nCount_Spatial_Whole.png", width = 10, height = 6)
ggsave("Vlnplot_nCount_Spatial_Whole.pdf", width = 10, height = 6)
VlnPlot(Whole_clean.merge, features = "nFeature_Spatial", group.by = "GEMgroup", pt.size = 0) + NoLegend()
ggsave("Vlnplot_nFeature_Spatial_Whole.png", width = 10, height = 6)
ggsave("Vlnplot_nFeature_Spatial_Whole.pdf", width = 10, height = 6)
VlnPlot(Whole_clean.merge, features = "percentMito", group.by = "GEMgroup", pt.size = 0) + NoLegend()
ggsave("Vlnplot_percentMito_Whole.png", width = 10, height = 6)
ggsave("Vlnplot_percentMito_Whole.pdf", width = 10, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_2") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_Y_B1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_Y_B1.pdf", width = 5, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_3") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_Y_C1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_Y_C1.pdf", width = 5, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_5") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_M_A1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_M_A1.pdf", width = 5, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_6") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_M_B1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_M_B1.pdf", width = 5, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_9") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_O_1_A1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_O_1_A1.pdf", width = 5, height = 6)

SpatialFeaturePlot(Whole_clean.merge, features = "nCount_Spatial", crop = F, 
                   images = "slice1_10") + 
  theme(legend.position = "right")
ggsave("SpatialFeaturePlot_nCount_Spatial_O_1_B1.png", width = 5, height = 6)
ggsave("SpatialFeaturePlot_nCount_Spatial_O_1_B1.pdf", width = 5, height = 6)

######################################################################################################
# SCT
######################################################################################################
DefaultAssay(Whole_clean.merge) <- "SCT"
VariableFeatures(Whole_clean.merge) <- c(VariableFeatures(Y_A1_clean), VariableFeatures(Y_B1_clean), 
                                         VariableFeatures(Y_C1_clean), VariableFeatures(Y_D1_clean), 
                                         VariableFeatures(M_A1_clean), VariableFeatures(M_B1_clean), 
                                         VariableFeatures(M_C1_clean), VariableFeatures(M_D1_clean), 
                                         VariableFeatures(O_1_A1_clean), VariableFeatures(O_1_B1_clean), 
                                         VariableFeatures(O_1_C1_clean), VariableFeatures(O_1_D1_clean), 
                                         VariableFeatures(O_2_A1_clean), VariableFeatures(O_2_B1_clean), 
                                         VariableFeatures(O_2_C1_clean), VariableFeatures(O_2_D1_clean))
Whole_clean.merge <- RunPCA(Whole_clean.merge, verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:50) %>% 
  RunUMAP(reduction = "pca", dims = 1:50) %>% 
  RunTSNE(reduction = "pca", dims = 1:50) %>% 
  FindClusters(verbose = FALSE)

DimPlot(Whole_clean.merge, reduction = "umap", label = T)
ggsave("UMAP_Visium_Whole_seurat_clusters.png", width = 7, height = 5.5)
ggsave("UMAP_Visium_Whole_seurat_clusters.pdf", width = 7, height = 5.5)

SpatialDimPlot(Whole_clean.merge, crop = F, images = "slice1_3")
ggsave("SpatialDimPlot_Whole_Y_C1.png", width = 5, height = 6)
ggsave("SpatialDimPlot_Whole_Y_C1.pdf", width = 5, height = 6)
SpatialDimPlot(Whole_clean.merge, crop = F, images = "slice1_5")
ggsave("SpatialDimPlot_Whole_M_A1.png", width = 5, height = 6)
ggsave("SpatialDimPlot_Whole_M_A1.pdf", width = 5, height = 6)
SpatialDimPlot(Whole_clean.merge, crop = F, images = "slice1_9")
ggsave("SpatialDimPlot_Whole_O_1_A1.png", width = 5, height = 6)
ggsave("SpatialDimPlot_Whole_O_1_A1.pdf", width = 5, height = 6)

######################################################################################################
# Annotation
######################################################################################################
Whole_clean.merge$Anatomy <- recode(Whole_clean.merge$seurat_clusters, 
                                    "0" = "Midbrain", 
                                    "1" = "Thalamus", 
                                    "2" = "Cortex", 
                                    "3" = "White_matter", 
                                    "4" = "Hippocampus", 
                                    "5" = "Cortex", 
                                    "6" = "Cortex", 
                                    "7" = "Cortex", 
                                    "8" = "Hippocampus", 
                                    "9" = "Hippocampus", 
                                    "10" = "Cortex", 
                                    "11" = "Hypothalamus", 
                                    "12" = "Cortex",
                                    "13" = "Hippocampus", 
                                    "14" = "Hippocampus", 
                                    "15" = "White_matter", 
                                    "16" = "White_matter", 
                                    "17" = "Hippocampus", 
                                    "18" = "White_matter", 
                                    "19" = "Cortex", 
                                    "20" = "White_matter", 
                                    "21" = "Cortex", 
                                    "22" = "Midbrain", 
                                    "23" = "Midbrain", 
                                    "24" = "Hypothalamus", 
                                    "25"  = "White_matter", 
                                    "26" = "Midbrain")
# Adjust the level of metadata
my_levels <- c("Cortex", "Hippocampus", "Thalamus", "Hypothalamus", "Midbrain", "White_matter")
Whole_clean.merge$Anatomy <- factor(x = Whole_clean.merge$Anatomy, levels = my_levels)
Idents(Whole_clean.merge) <- "Anatomy"
DimPlot(Whole_clean.merge, reduction = "umap", 
        cols = c("#807dba", "#f16913", "#1d91c0", "#41ab5d", "#feb24c", "#ef3b2c"))
ggsave("UMAP_Visium_Whole.png", width = 7, height = 5.5)
ggsave("UMAP_Visium_Whole.pdf", width = 7, height = 5.5)

######################################################################################################
# Save results
######################################################################################################
saveRDS(Whole_clean.merge, "Whole_clean.merge.rds")


######################################################################################################
# Transfer the label of Hippocampus
######################################################################################################
## Young
matadata_Y_A1 <- read.csv("./Data/Visium coordinates/Y_A1_Hippocampus.csv")
names <- make.unique(matadata_Y_A1$Barcode)
rownames(matadata_Y_A1) <- names
colnames(matadata_Y_A1)[2] <- "Region"
Y_A1 <- AddMetaData(Y_A1, metadata = matadata_Y_A1)
Y_A1$Barcode <- NULL
# Subset
Idents(Y_A1) <- "Region"
Y_A1_Hippo <- subset(Y_A1, idents = c("ML", "GCL", "Hilus", "CA1"))

matadata_Y_B1 <- read.csv("./Data/Visium coordinates/Y_B1_Hippocampus.csv")
names <- make.unique(matadata_Y_B1$Barcode)
rownames(matadata_Y_B1) <- names
colnames(matadata_Y_B1)[2] <- "Region"
Y_B1 <- AddMetaData(Y_B1, metadata = matadata_Y_B1)
Y_B1$Barcode <- NULL
# Subset
Idents(Y_B1) <- "Region"
Y_B1_Hippo <- subset(Y_B1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                      "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

matadata_Y_C1 <- read.csv("./Data/Visium coordinates/Y_C1_Hippocampus.csv")
names <- make.unique(matadata_Y_C1$Barcode)
rownames(matadata_Y_C1) <- names
colnames(matadata_Y_C1)[2] <- "Region"
Y_C1 <- AddMetaData(Y_C1, metadata = matadata_Y_C1)
Y_C1$Barcode <- NULL
# Subset
Idents(Y_C1) <- "Region"
Y_C1_Hippo <- subset(Y_C1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1", "CA3", 
                                      "ML_2", "GCL_2", "Hilus_2"))

matadata_Y_D1 <- read.csv("./Data/Visium coordinates/Y_D1_Hippocampus.csv")
names <- make.unique(matadata_Y_D1$Barcode)
rownames(matadata_Y_D1) <- names
colnames(matadata_Y_D1)[2] <- "Region"
Y_D1 <- AddMetaData(Y_D1, metadata = matadata_Y_D1)
Y_D1$Barcode <- NULL
# Subset
Idents(Y_D1) <- "Region"
Y_D1_Hippo <- subset(Y_D1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                      "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

## Middle-age
matadata_M_A1 <- read.csv("./Data/Visium coordinates/M_A1_Hippocampus.csv")
names <- make.unique(matadata_M_A1$Barcode)
rownames(matadata_M_A1) <- names
colnames(matadata_M_A1)[2] <- "Region"
M_A1 <- AddMetaData(M_A1, metadata = matadata_M_A1)
M_A1$Barcode <- NULL
# Subset
Idents(M_A1) <- "Region"
M_A1_Hippo <- subset(M_A1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                      "ML_2", "GCL_2", "Hilus_2", 
                                      "CA1", "CA3"))

matadata_M_B1 <- read.csv("./Data/Visium coordinates/M_B1_Hippocampus.csv")
names <- make.unique(matadata_M_B1$Barcode)
rownames(matadata_M_B1) <- names
colnames(matadata_M_B1)[2] <- "Region"
M_B1 <- AddMetaData(M_B1, metadata = matadata_M_B1)
M_B1$Barcode <- NULL
# Subset
Idents(M_B1) <- "Region"
M_B1_Hippo <- subset(M_B1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                      "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

matadata_M_C1 <- read.csv("./Data/Visium coordinates/M_C1_Hippocampus.csv")
names <- make.unique(matadata_M_C1$Barcode)
rownames(matadata_M_C1) <- names
colnames(matadata_M_C1)[2] <- "Region"
M_C1 <- AddMetaData(M_C1, metadata = matadata_M_C1)
M_C1$Barcode <- NULL
# Subset
Idents(M_C1) <- "Region"
M_C1_Hippo <- subset(M_C1, idents = c("ML", "GCL", "Hilus", "CA1", "CA3"))

matadata_M_D1 <- read.csv("./Data/Visium coordinates/M_D1_Hippocampus.csv")
names <- make.unique(matadata_M_D1$Barcode)
rownames(matadata_M_D1) <- names
colnames(matadata_M_D1)[2] <- "Region"
M_D1 <- AddMetaData(M_D1, metadata = matadata_M_D1)
M_D1$Barcode <- NULL
# Subset
Idents(M_D1) <- "Region"
M_D1_Hippo <- subset(M_D1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                      "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

## Old
matadata_O_1_A1 <- read.csv("./Data/Visium coordinates/O_1_A1_Hippocampus.csv")
names <- make.unique(matadata_O_1_A1$Barcode)
rownames(matadata_O_1_A1) <- names
colnames(matadata_O_1_A1)[2] <- "Region"
O_1_A1 <- AddMetaData(O_1_A1, metadata = matadata_O_1_A1)
O_1_A1$Barcode <- NULL
# Subset
Idents(O_1_A1) <- "Region"
O_1_A1_Hippo <- subset(O_1_A1, idents = c("ML", "GCL", "Hilus", "CA1", "CA3"))

matadata_O_1_B1 <- read.csv("./Data/Visium coordinates/O_1_B1_Hippocampus.csv")
names <- make.unique(matadata_O_1_B1$Barcode)
rownames(matadata_O_1_B1) <- names
colnames(matadata_O_1_B1)[2] <- "Region"
O_1_B1 <- AddMetaData(O_1_B1, metadata = matadata_O_1_B1)
O_1_B1$Barcode <- NULL
# Subset
Idents(O_1_B1) <- "Region"
O_1_B1_Hippo <- subset(O_1_B1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                          "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

matadata_O_1_C1 <- read.csv("./Data/Visium coordinates/O_1_C1_Hippocampus.csv")
names <- make.unique(matadata_O_1_C1$Barcode)
rownames(matadata_O_1_C1) <- names
colnames(matadata_O_1_C1)[2] <- "Region"
O_1_C1 <- AddMetaData(O_1_C1, metadata = matadata_O_1_C1)
O_1_C1$Barcode <- NULL
# Subset
Idents(O_1_C1) <- "Region"
O_1_C1_Hippo <- subset(O_1_C1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1", "CA2", "CA3", 
                                          "ML_2", "GCL_2", "Hilus_2"))

matadata_O_1_D1 <- read.csv("./Data/Visium coordinates/O_1_D1_Hippocampus.csv")
names <- make.unique(matadata_O_1_D1$Barcode)
rownames(matadata_O_1_D1) <- names
colnames(matadata_O_1_D1)[2] <- "Region"
O_1_D1 <- AddMetaData(O_1_D1, metadata = matadata_O_1_D1)
O_1_D1$Barcode <- NULL
# Subset
Idents(O_1_D1) <- "Region"
O_1_D1_Hippo <- subset(O_1_D1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                          "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

matadata_O_2_A1 <- read.csv("./Data/Visium coordinates/O_2_A1_Hippocampus.csv")
names <- make.unique(matadata_O_2_A1$Barcode)
rownames(matadata_O_2_A1) <- names
colnames(matadata_O_2_A1)[2] <- "Region"
O_2_A1 <- AddMetaData(O_2_A1, metadata = matadata_O_2_A1)
O_2_A1$Barcode <- NULL
# Subset
Idents(O_2_A1) <- "Region"
O_2_A1_Hippo <- subset(O_2_A1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1", "CA3", 
                                          "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_B1 <- read.csv("./Data/Visium coordinates/O_2_B1_Hippocampus.csv")
names <- make.unique(matadata_O_2_B1$Barcode)
rownames(matadata_O_2_B1) <- names
colnames(matadata_O_2_B1)[2] <- "Region"
O_2_B1 <- AddMetaData(O_2_B1, metadata = matadata_O_2_B1)
O_2_B1$Barcode <- NULL
# Subset
Idents(O_2_B1) <- "Region"
O_2_B1_Hippo <- subset(O_2_B1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                          "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))

matadata_O_2_C1 <- read.csv("./Data/Visium coordinates/O_2_C1_Hippocampus.csv")
names <- make.unique(matadata_O_2_C1$Barcode)
rownames(matadata_O_2_C1) <- names
colnames(matadata_O_2_C1)[2] <- "Region"
O_2_C1 <- AddMetaData(O_2_C1, metadata = matadata_O_2_C1)
O_2_C1$Barcode <- NULL
# Subset
Idents(O_2_C1) <- "Region"
O_2_C1_Hippo <- subset(O_2_C1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1", "CA2", "CA3", 
                                          "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_D1 <- read.csv("./Data/Visium coordinates/O_2_D1_Hippocampus.csv")
names <- make.unique(matadata_O_2_D1$Barcode)
rownames(matadata_O_2_D1) <- names
colnames(matadata_O_2_D1)[2] <- "Region"
O_2_D1 <- AddMetaData(O_2_D1, metadata = matadata_O_2_D1)
O_2_D1$Barcode <- NULL
# Subset
Idents(O_2_D1) <- "Region"
O_2_D1_Hippo <- subset(O_2_D1, idents = c("ML_1", "GCL_1", "Hilus_1", "CA1_1", "CA2_1", "CA3_1", 
                                          "ML_2", "GCL_2", "Hilus_2", "CA1_2", "CA2_2", "CA3_2"))
######################################################################################################
# Transfer the label of DG
######################################################################################################
## Young
matadata_Y_A1 <- read.csv("./Data/Visium coordinates/Y_A1_Hippocampus.csv")
names <- make.unique(matadata_Y_A1$Barcode)
rownames(matadata_Y_A1) <- names
colnames(matadata_Y_A1)[2] <- "Region"
Y_A1 <- AddMetaData(Y_A1, metadata = matadata_Y_A1)
Y_A1$Barcode <- NULL
# Subset
Idents(Y_A1) <- "Region"
Y_A1_DG <- subset(Y_A1, idents = c("ML", "GCL", "Hilus"))

matadata_Y_B1 <- read.csv("./Data/Visium coordinates/Y_B1_Hippocampus.csv")
names <- make.unique(matadata_Y_B1$Barcode)
rownames(matadata_Y_B1) <- names
colnames(matadata_Y_B1)[2] <- "Region"
Y_B1 <- AddMetaData(Y_B1, metadata = matadata_Y_B1)
Y_B1$Barcode <- NULL
# Subset
Idents(Y_B1) <- "Region"
Y_B1_DG <- subset(Y_B1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

matadata_Y_C1 <- read.csv("./Data/Visium coordinates/Y_C1_Hippocampus.csv")
names <- make.unique(matadata_Y_C1$Barcode)
rownames(matadata_Y_C1) <- names
colnames(matadata_Y_C1)[2] <- "Region"
Y_C1 <- AddMetaData(Y_C1, metadata = matadata_Y_C1)
Y_C1$Barcode <- NULL
# Subset
Idents(Y_C1) <- "Region"
Y_C1_DG <- subset(Y_C1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

matadata_Y_D1 <- read.csv("./Data/Visium coordinates/Y_D1_Hippocampus.csv")
names <- make.unique(matadata_Y_D1$Barcode)
rownames(matadata_Y_D1) <- names
colnames(matadata_Y_D1)[2] <- "Region"
Y_D1 <- AddMetaData(Y_D1, metadata = matadata_Y_D1)
Y_D1$Barcode <- NULL
# Subset
Idents(Y_D1) <- "Region"
Y_D1_DG <- subset(Y_D1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

## Middle-age
matadata_M_A1 <- read.csv("./Data/Visium coordinates/M_A1_Hippocampus.csv")
names <- make.unique(matadata_M_A1$Barcode)
rownames(matadata_M_A1) <- names
colnames(matadata_M_A1)[2] <- "Region"
M_A1 <- AddMetaData(M_A1, metadata = matadata_M_A1)
M_A1$Barcode <- NULL
# Subset
Idents(M_A1) <- "Region"
M_A1_DG <- subset(M_A1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

matadata_M_B1 <- read.csv("./Data/Visium coordinates/M_B1_Hippocampus.csv")
names <- make.unique(matadata_M_B1$Barcode)
rownames(matadata_M_B1) <- names
colnames(matadata_M_B1)[2] <- "Region"
M_B1 <- AddMetaData(M_B1, metadata = matadata_M_B1)
M_B1$Barcode <- NULL
# Subset
Idents(M_B1) <- "Region"
M_B1_DG <- subset(M_B1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

matadata_M_C1 <- read.csv("./Data/Visium coordinates/M_C1_Hippocampus.csv")
names <- make.unique(matadata_M_C1$Barcode)
rownames(matadata_M_C1) <- names
colnames(matadata_M_C1)[2] <- "Region"
M_C1 <- AddMetaData(M_C1, metadata = matadata_M_C1)
M_C1$Barcode <- NULL
# Subset
Idents(M_C1) <- "Region"
M_C1_DG <- subset(M_C1, idents = c("ML", "GCL", "Hilus"))

matadata_M_D1 <- read.csv("./Data/Visium coordinates/M_D1_Hippocampus.csv")
names <- make.unique(matadata_M_D1$Barcode)
rownames(matadata_M_D1) <- names
colnames(matadata_M_D1)[2] <- "Region"
M_D1 <- AddMetaData(M_D1, metadata = matadata_M_D1)
M_D1$Barcode <- NULL
# Subset
Idents(M_D1) <- "Region"
M_D1_DG <- subset(M_D1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                   "ML_2", "GCL_2", "Hilus_2"))

## Old
matadata_O_1_A1 <- read.csv("./Data/Visium coordinates/O_1_A1_Hippocampus.csv")
names <- make.unique(matadata_O_1_A1$Barcode)
rownames(matadata_O_1_A1) <- names
colnames(matadata_O_1_A1)[2] <- "Region"
O_1_A1 <- AddMetaData(O_1_A1, metadata = matadata_O_1_A1)
O_1_A1$Barcode <- NULL
# Subset
Idents(O_1_A1) <- "Region"
O_1_A1_DG <- subset(O_1_A1, idents = c("ML", "GCL", "Hilus"))

matadata_O_1_B1 <- read.csv("./Data/Visium coordinates/O_1_B1_Hippocampus.csv")
names <- make.unique(matadata_O_1_B1$Barcode)
rownames(matadata_O_1_B1) <- names
colnames(matadata_O_1_B1)[2] <- "Region"
O_1_B1 <- AddMetaData(O_1_B1, metadata = matadata_O_1_B1)
O_1_B1$Barcode <- NULL
# Subset
Idents(O_1_B1) <- "Region"
O_1_B1_DG <- subset(O_1_B1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_1_C1 <- read.csv("./Data/Visium coordinates/O_1_C1_Hippocampus.csv")
names <- make.unique(matadata_O_1_C1$Barcode)
rownames(matadata_O_1_C1) <- names
colnames(matadata_O_1_C1)[2] <- "Region"
O_1_C1 <- AddMetaData(O_1_C1, metadata = matadata_O_1_C1)
O_1_C1$Barcode <- NULL
# Subset
Idents(O_1_C1) <- "Region"
O_1_C1_DG <- subset(O_1_C1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_1_D1 <- read.csv("./Data/Visium coordinates/O_1_D1_Hippocampus.csv")
names <- make.unique(matadata_O_1_D1$Barcode)
rownames(matadata_O_1_D1) <- names
colnames(matadata_O_1_D1)[2] <- "Region"
O_1_D1 <- AddMetaData(O_1_D1, metadata = matadata_O_1_D1)
O_1_D1$Barcode <- NULL
# Subset
Idents(O_1_D1) <- "Region"
O_1_D1_DG <- subset(O_1_D1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_A1 <- read.csv("./Data/Visium coordinates/O_2_A1_Hippocampus.csv")
names <- make.unique(matadata_O_2_A1$Barcode)
rownames(matadata_O_2_A1) <- names
colnames(matadata_O_2_A1)[2] <- "Region"
O_2_A1 <- AddMetaData(O_2_A1, metadata = matadata_O_2_A1)
O_2_A1$Barcode <- NULL
# Subset
Idents(O_2_A1) <- "Region"
O_2_A1_DG <- subset(O_2_A1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_B1 <- read.csv("./Data/Visium coordinates/O_2_B1_Hippocampus.csv")
names <- make.unique(matadata_O_2_B1$Barcode)
rownames(matadata_O_2_B1) <- names
colnames(matadata_O_2_B1)[2] <- "Region"
O_2_B1 <- AddMetaData(O_2_B1, metadata = matadata_O_2_B1)
O_2_B1$Barcode <- NULL
# Subset
Idents(O_2_B1) <- "Region"
O_2_B1_DG <- subset(O_2_B1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_C1 <- read.csv("./Data/Visium coordinates/O_2_C1_Hippocampus.csv")
names <- make.unique(matadata_O_2_C1$Barcode)
rownames(matadata_O_2_C1) <- names
colnames(matadata_O_2_C1)[2] <- "Region"
O_2_C1 <- AddMetaData(O_2_C1, metadata = matadata_O_2_C1)
O_2_C1$Barcode <- NULL
# Subset
Idents(O_2_C1) <- "Region"
O_2_C1_DG <- subset(O_2_C1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))

matadata_O_2_D1 <- read.csv("./Data/Visium coordinates/O_2_D1_Hippocampus.csv")
names <- make.unique(matadata_O_2_D1$Barcode)
rownames(matadata_O_2_D1) <- names
colnames(matadata_O_2_D1)[2] <- "Region"
O_2_D1 <- AddMetaData(O_2_D1, metadata = matadata_O_2_D1)
O_2_D1$Barcode <- NULL
# Subset
Idents(O_2_D1) <- "Region"
O_2_D1_DG <- subset(O_2_D1, idents = c("ML_1", "GCL_1", "Hilus_1", 
                                       "ML_2", "GCL_2", "Hilus_2"))








##== Integration with Allen Brain Atlas HP ==##
Hippo_All <- readRDS("Data/SeuratObject/Visium/Processed/Hippo_All.rds")
Hippo_All <- SCTransform(Hippo_All, assay = "Spatial", verbose = F) %>% 
  RunPCA(verbose = F) ## Default: ncells = 5000, variable.features.n = 3000

Hippo <- readRDS("Data/SeuratObject/Chromium/Integration/AllenBrainAtlas_Hippo.rds")
Hippo <- SCTransform(Hippo, verbose = F) %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30)
Hippo$Celltype <- recode(Hippo$subclass_label, 
                         "Astro" = "Astrocyte", 
                         "Oligo" = "Oligodendrocyte", 
                         "Micro-PVM" = "Micro-PVM", 
                         "CA1-ProS" = "CA1", 
                         "CA2-IG-FC" = "CA2", 
                         "CA3" = "CA3", 
                         "DG" = "DG", 
                         "Lamp5" = "Interneuron", 
                         "Pvalb" = "Interneuron", 
                         "Sncg" = "Interneuron", 
                         "Sst" = "Interneuron", 
                         "Sst Chodl" = "Interneuron", 
                         "Vip" = "Interneuron", 
                         "L2/3 IT ENTl" = "RHP", 
                         "L2/3 IT RHP" = "RHP", 
                         "L6 CT CTX" = "RHP", 
                         "L6b CTX" = "RHP", 
                         "NP SUB"  = "RHP", 
                         "SUB-ProS" = "RHP")
Idents(Hippo) <- "Celltype"
DimPlot(Hippo, reduction = "umap", 
        cols = c("#dd3497", "#993404", "#fe9929", "#ef3b2c", "#ec7014", "#fb6a4a", 
                 "#969696", "#a8ddb5", "#4eb3d3", "#1d91c0", "#9e9ac8"))
######################################################################################################
# Identify anchors using the reference data
######################################################################################################
anchors <- FindTransferAnchors(reference = Hippo, query = Hippo_All, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Hippo$Celltype, prediction.assay = TRUE, 
                                  weight.reduction = Hippo_All[["pca"]], dims = 1:30, k.weight = 30)
Hippo_All[["predictions"]] <- predictions.assay
