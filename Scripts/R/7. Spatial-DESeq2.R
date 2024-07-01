library(DESeq2)

##== Data preparation ==##
######################################################################################################
# Read SeuratObject
######################################################################################################
# Bring in Seurat object
seurat <- readRDS("./Data/SeuratObject/Visium/Processed/Hippo_All.rds")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$Spatial@counts
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat$Region)
metadata$group_id <- factor(seurat$Age)
metadata$sample_id <- factor(seurat$GEMgroup)

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
save.image(file = 'DESeq2_Visium_Hippo_All.RData')






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
# Subsetting DG
######################################################################################################
# Generate vector of cluster IDs
clusters <- unique(metadata$cluster_id)
clusters
clusters[1]

# Subset the metadata to only the DG
cluster_metadata <- metadata[which(metadata$cluster_id == "DG"), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the DG
counts <- pb[["DG"]]
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
  scale_color_manual(breaks = c("Young", "Middle", "Old"), 
                     values = c("#c6dbef", "#4292c6", "#08306b")) + 
  theme_bw()
ggsave("DESeq2_PCA_Visium_DG_All.png", width = 5, height = 2.5)
ggsave("DESeq2_PCA_Visium_DG_All.pdf", width = 5, height = 2.5)

######################################################################################################
# Hierarchical clustering
######################################################################################################
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
pheatmap(rld_cor, 
         annotation = cluster_metadata[, c("group_id"), drop = F], 
         color = colors, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         filename = "DESeq2_Visium_Heatmap_DG.pdf",
         width = 6.2, height = 5)

######################################################################################################
# Run DESeq2 differential expression analysis
######################################################################################################
dds <- DESeq(dds)
# Plot dispersion estimates
pdf("DispersionPlot_Visium_DG_All.pdf", width = 5, height = 4)
plotDispEsts(dds)
dev.off()




##== Visualization ==##
######################################################################################################
# VolcanoPlot 
######################################################################################################
# OY
df <- results(dds, contrast = c('group_id', 'Old', 'Young'))
head(df, 3)
df_selected <- as.data.frame(df) %>% 
  drop_na()
df_selected$gene <- rownames(df_selected)
df_selected$Change <- ifelse(df_selected$padj < 0.05,
                             ifelse(abs(df_selected$log2FoldChange) > 0,
                                    ifelse(df_selected$log2FoldChange < 0, 'Down', 'Up'), 'Stable'), 'Stable')
## Check
top50 <- df_selected[order(-log10(df_selected$padj), decreasing = T), ]
top50 <- top50 %>% 
  filter(Change != "Stable") 
top50 <- top50[1:50, ]
## Volcano plot
ggplot(df_selected, aes(x = log2FoldChange, y = -log10(padj), fill = "Change")) +
  geom_point(aes(colour = Change)) + 
  scale_color_manual(values = c("#6baed6", "lightgrey", "firebrick")) + 
  ggtitle("DG Old vs. Young") +
  xlab(log[2]~FC) + 
  ylab(-log[10]~p_adj_val) + 
  geom_hline(aes(yintercept = -log10(0.05)), color  = "darkgrey", linetype = "dashed") +
  geom_vline(aes(xintercept = log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_vline(aes(xintercept = -log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_text_repel(data = head(top50, 50), aes(label = gene)) + 
  theme_classic() + 
  NoLegend()
ggsave("DESeq2_Visium_Volcanoplot_OY_DG.png", width = 9, height = 8)
ggsave("DESeq2_Visium_Volcanoplot_OY_DG.pdf", width = 9, height = 8)

# MY
df <- results(dds, contrast = c('group_id', 'Middle', 'Young'))
head(df, 3)
df_selected <- as.data.frame(df) %>% 
  drop_na()
df_selected$gene <- rownames(df_selected)
df_selected$Change <- ifelse(df_selected$padj < 0.05,
                             ifelse(abs(df_selected$log2FoldChange) > 0,
                                    ifelse(df_selected$log2FoldChange < 0, 'Down', 'Up'), 'Stable'), 'Stable')
## Check
top50 <- df_selected[order(-log10(df_selected$padj), decreasing = T), ]
top50 <- top50 %>% 
  filter(Change != "Stable") 
top50 <- top50[1:50, ]
## Volcano plot
ggplot(df_selected, aes(x = log2FoldChange, y = -log10(padj), fill = "Change")) +
  geom_point(aes(colour = Change)) + 
  scale_color_manual(values = c("#6baed6", "lightgrey", "firebrick")) + 
  ggtitle("DG Middle vs. Young") +
  xlab(log[2]~FC) + 
  ylab(-log[10]~p_adj_val) + 
  geom_hline(aes(yintercept = -log10(0.05)), color  = "darkgrey", linetype = "dashed") +
  geom_vline(aes(xintercept = log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_vline(aes(xintercept = -log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_text_repel(data = head(top50, 50), aes(label = gene)) + 
  theme_classic() + 
  NoLegend()
ggsave("DESeq2_Visium_Volcanoplot_MY_DG.png", width = 9, height = 8)
ggsave("DESeq2_Visium_Volcanoplot_MY_DG.pdf", width = 9, height = 8)

# OM
df <- results(dds, contrast = c('group_id', 'Old', 'Middle'))
head(df, 3)
df_selected <- as.data.frame(df) %>% 
  drop_na()
df_selected$gene <- rownames(df_selected)
df_selected$Change <- ifelse(df_selected$padj < 0.05,
                             ifelse(abs(df_selected$log2FoldChange) > 0,
                                    ifelse(df_selected$log2FoldChange < 0, 'Down', 'Up'), 'Stable'), 'Stable')
## Check
top50 <- df_selected[order(-log10(df_selected$padj), decreasing = T), ]
top50 <- top50 %>% 
  filter(Change != "Stable") 
top50 <- top50[1:50, ]
## Volcano plot
ggplot(df_selected, aes(x = log2FoldChange, y = -log10(padj), fill = "Change")) +
  geom_point(aes(colour = Change)) + 
  scale_color_manual(values = c("lightgrey", "firebrick")) + 
  ggtitle("DG Old vs. Middle") +
  xlab(log[2]~FC) + 
  ylab(-log[10]~p_adj_val) + 
  geom_hline(aes(yintercept = -log10(0.05)), color  = "darkgrey", linetype = "dashed") +
  geom_vline(aes(xintercept = log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_vline(aes(xintercept = -log2(1.5)), color  = "darkgrey", linetype = "dashed") + 
  geom_text_repel(data = head(top50, 50), aes(label = gene)) + 
  theme_classic() + 
  NoLegend()
ggsave("DESeq2_Visium_Volcanoplot_OM_DG.png", width = 9, height = 8)
ggsave("DESeq2_Visium_Volcanoplot_OM_DG.pdf", width = 9, height = 8)

######################################################################################################
# Heatmap
######################################################################################################
# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- df %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < 0.05) %>%
  dplyr::arrange(padj)

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- data.frame(normalized_counts)
colnames(normalized_counts) <- gsub("X", "", colnames(normalized_counts))
sig_norm <- normalized_counts %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

heat_colors <- brewer.pal(6, "RdPu")
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("group_id", "cluster_id")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         filename = "DESeq2_Visium_Heatmap_DG_DEG.pdf",
         width = 6.2, height = 5)  





