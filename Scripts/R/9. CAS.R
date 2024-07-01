# Load library
rm(list=ls())
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(UpSetR)
library(RColorBrewer)


##== Up-regulation ==##
######################################################################################################
# Load data
######################################################################################################
# Whole DG
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Visium_DG_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "6")
df_2 <- df %>% 
  filter(cluster == "8")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
DG_Up <- c(names_1, names_2)
# qNSC
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_qNSC_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "2")
df_2 <- df %>% 
  filter(cluster == "4")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
qNSC_Up <- c(names_1, names_2)
# Astrocyte
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Astrocyte_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "4")
df_2 <- df %>% 
  filter(cluster == "7")
df_3 <- df %>% 
  filter(cluster == "8")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
Astrocyte_Up <- c(names_1, names_2, names_3)
# Microglia
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Microglia_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "1")
df_2 <- df %>% 
  filter(cluster == "4")
df_3 <- df %>% 
  filter(cluster == "7")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
Microglia_Up <- c(names_1, names_2, names_3)
# Endothelial
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Endothelial_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "5")
df_2 <- df %>% 
  filter(cluster == "6")
df_3 <- df %>% 
  filter(cluster == "8")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
EC_Up <- c(names_1, names_2, names_3)

######################################################################################################
# UpsetR plot
######################################################################################################
set_list <- list("Astrocyte" = Astrocyte_Up, 
                 "qNSC" = qNSC_Up, 
                 "Endothelial" = EC_Up, 
                 "Microglia" = Microglia_Up, 
                 "DG" = DG_Up)
# Plotting
upset(fromList(set_list), 
      order.by = c("freq"), 
      sets = c("Astrocyte", "qNSC", "Endothelial", "Microglia", "DG"), 
      sets.bar.color = brewer.pal(5, "Set3"), 
      matrix.color = brewer.pal(4, "Set1")[2])
# Save the figure
p1 <- upset(fromList(set_list), 
            order.by = c("freq"), 
            sets = c("Astrocyte", "qNSC", "Endothelial", "Microglia", "DG"), 
            sets.bar.color = brewer.pal(5, "Set3"), 
            matrix.color = brewer.pal(4, "Set1")[2])
pdf("UpSet_GenericAgingSignature_Up.pdf", width = 9, height = 4.5, useDingbats = FALSE)
p1
dev.off()

######################################################################################################
# Gene list
######################################################################################################
Common_Up_1 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_1 <- Common_Up_1[Common_Up_1 %in% Astrocyte_Up]

Common_Up_2 <- DG_Up[DG_Up %in% Astrocyte_Up]
Common_Up_2 <- Common_Up_2[Common_Up_2 %in% Microglia_Up]

Common_Up_3 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_3 <- Common_Up_3[Common_Up_3 %in% Microglia_Up]

Common_Up_4 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_4 <- Common_Up_4[Common_Up_4 %in% Microglia_Up]
Common_Up_4 <- Common_Up_4[Common_Up_4 %in% Astrocyte_Up]

Common_Up_5 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_5 <- Common_Up_5[Common_Up_5 %in% qNSC_Up]

Common_Up_6 <- DG_Up[DG_Up %in% Astrocyte_Up]
Common_Up_6 <- Common_Up_6[Common_Up_6 %in% qNSC_Up]

Common_Up_7 <- DG_Up[DG_Up %in% Microglia_Up]
Common_Up_7 <- Common_Up_7[Common_Up_7 %in% qNSC_Up]

Common_Up_8 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_8 <- Common_Up_8[Common_Up_8 %in% qNSC_Up]
Common_Up_8 <- Common_Up_8[Common_Up_8 %in% Astrocyte_Up]

Common_Up_9 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_9 <- Common_Up_9[Common_Up_9 %in% qNSC_Up]
Common_Up_9 <- Common_Up_9[Common_Up_9 %in% Microglia_Up]

Common_Up_10 <- DG_Up[DG_Up %in% EC_Up]
Common_Up_10 <- Common_Up_10[Common_Up_10 %in% qNSC_Up]
Common_Up_10 <- Common_Up_10[Common_Up_10 %in% Microglia_Up]
Common_Up_10 <- Common_Up_10[Common_Up_10 %in% Astrocyte_Up]

Common_Up_11 <- DG_Up[DG_Up %in% Astrocyte_Up]
Common_Up_11 <- Common_Up_11[Common_Up_11 %in% qNSC_Up]
Common_Up_11 <- Common_Up_11[Common_Up_11 %in% Microglia_Up]

Common_Up <- unique(c(Common_Up_1, Common_Up_2, Common_Up_3, Common_Up_4, Common_Up_5, Common_Up_6, Common_Up_7, 
                      Common_Up_8, Common_Up_9, Common_Up_10, Common_Up_11))
# Save data
write.csv(Common_Up, "GenericAgingSignature_Up.csv") # Used for GSEA


##== Down-regulation ==##
######################################################################################################
# Load data
######################################################################################################
# Whole DG
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Visium_DG_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "4")
df_2 <- df %>% 
  filter(cluster == "5")
df_3 <- df %>% 
  filter(cluster == "7")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
DG_Down <- c(names_1, names_2, names_3)
# qNSC
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_qNSC_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "5")
df_2 <- df %>% 
  filter(cluster == "6")
df_3 <- df %>% 
  filter(cluster == "8")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
qNSC_Down <- c(names_1, names_2, names_3)
# Astrocyte
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Astrocyte_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "2")
df_2 <- df %>% 
  filter(cluster == "5")
df_3 <- df %>% 
  filter(cluster == "6")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
Astrocyte_Down <- c(names_1, names_2, names_3)
# Microglia
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Microglia_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "2")
df_2 <- df %>% 
  filter(cluster == "6")
df_3 <- df %>% 
  filter(cluster == "8")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
Microglia_Down <- c(names_1, names_2, names_3)
# Endothelial
df <- read.csv("./Data/DEG/Mfuzz/GeneList/Mfuzz_Endothelial_HVG.csv")
df <- df[, -c(1, 3:10)]
df_1 <- df %>% 
  filter(cluster == "1")
df_2 <- df %>% 
  filter(cluster == "4")
df_3 <- df %>% 
  filter(cluster == "7")
names_1 <- make.unique(df_1$gene)
names_2 <- make.unique(df_2$gene)
names_3 <- make.unique(df_3$gene)
EC_Down <- c(names_1, names_2, names_3)

######################################################################################################
# UpsetR plot
######################################################################################################
set_list <- list("Astrocyte" = Astrocyte_Down, 
                 "qNSC" = qNSC_Down, 
                 "Endothelial" = EC_Down, 
                 "Microglia" = Microglia_Down, 
                 "DG" = DG_Down)
# Plotting
upset(fromList(set_list), 
      order.by = c("freq"), 
      sets = c("Astrocyte", "qNSC", "Endothelial", "Microglia", "DG"), 
      sets.bar.color = brewer.pal(5, "Set3"), 
      matrix.color = brewer.pal(4, "Set1")[2])
# Save the figure
p2 <- upset(fromList(set_list), 
            order.by = c("freq"), 
            sets = c("Astrocyte", "qNSC", "Endothelial", "Microglia", "DG"), 
            sets.bar.color = brewer.pal(5, "Set3"), 
            matrix.color = brewer.pal(4, "Set1")[2])
pdf("UpSet_GenericAgingSignature_Down.pdf", width = 9, height = 4.5, useDingbats = FALSE)
p2
dev.off()

######################################################################################################
# Gene list
######################################################################################################
Common_Down_1 <- DG_Down[DG_Down %in% Microglia_Down]
Common_Down_1 <- Common_Down_1[Common_Down_1 %in% Astrocyte_Down]

Common_Down_2 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_2 <- Common_Down_2[Common_Down_2 %in% Astrocyte_Down]

Common_Down_3 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_3 <- Common_Down_3[Common_Down_3 %in% Microglia_Down]

Common_Down_4 <- DG_Down[DG_Down %in% EC_Down]
Common_Down_4 <- Common_Down_4[Common_Down_4 %in% Astrocyte_Down]

Common_Down_5 <- DG_Down[DG_Down %in% EC_Down]
Common_Down_5 <- Common_Down_5[Common_Down_5 %in% Astrocyte_Down]
Common_Down_5 <- Common_Down_5[Common_Down_5 %in% Microglia_Down]

Common_Down_6 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_6 <- Common_Down_6[Common_Down_6 %in% EC_Down]

Common_Down_7 <- DG_Down[DG_Down %in% Microglia_Down]
Common_Down_7 <- Common_Down_7[Common_Down_7 %in% EC_Down]

Common_Down_8 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_8 <- Common_Down_8[Common_Down_8 %in% Astrocyte_Down]
Common_Down_8 <- Common_Down_8[Common_Down_8 %in% Microglia_Down]

Common_Down_9 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_9 <- Common_Down_9[Common_Down_9 %in% Astrocyte_Down]
Common_Down_9 <- Common_Down_9[Common_Down_9 %in% EC_Down]

Common_Down_10 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_10 <- Common_Down_10[Common_Down_10 %in% Astrocyte_Down]
Common_Down_10 <- Common_Down_10[Common_Down_10 %in% Microglia_Down]
Common_Down_10 <- Common_Down_10[Common_Down_10 %in% EC_Down]

Common_Down_11 <- DG_Down[DG_Down %in% qNSC_Down]
Common_Down_11 <- Common_Down_11[Common_Down_11 %in% Microglia_Down]
Common_Down_11 <- Common_Down_11[Common_Down_11 %in% EC_Down]

Common_Down <- unique(c(Common_Down_1, Common_Down_2, Common_Down_3, Common_Down_4, Common_Down_5, Common_Down_6, 
                 Common_Down_7, Common_Down_8, Common_Down_9, Common_Down_10, Common_Down_11))
# Save data
write.csv(Common_Down, "CommonAgingSignature_Down.csv") # Used for GSEA


