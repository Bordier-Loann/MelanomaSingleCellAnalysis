# Set the working directory to the specified path
setwd("~/Projet_Sante/Code")

# Load the raw data file containing single-cell expression data
counts <- read.table('GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz', header = TRUE, fill = TRUE, skip = 1)

# Display the first few rows of the expression data for inspection
head(counts)

# Install the 'readr' package if not already installed
install.packages('readr')

# Load the 'readr' package for reading and processing data
library(readr)

# Clean the raw metadata file containing sample information
metadata <- read_tsv('GSE120575_patient_ID_single_cells.txt.gz', skip = 19, col_types = cols(
  'Sample_name' = col_character(),
  'title' = col_character(),
  'source_name' = col_character(),
  'organism' = col_character(),
  'characteristics: patinet ID (Pre=baseline; Post= on treatment)' = col_character(),
  'characteristics: response' = col_character(),
  'characteristics: therapy' = col_character()
))

# Display the first few rows of the cleaned metadata for inspection
head(metadata)

# Remove columns with all NA values from metadata
metadata <- metadata[, colSums(!is.na(metadata)) > 0]

# Print the head of the metadata after cleaning
print(head(metadata))

# Retain only the first 16291 rows of the metadata
metadata <- metadata[1:16291, ]

# Display the tail of the metadata for inspection
tail(metadata)

# Load necessary libraries for downstream analysis
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

# Rename metadata columns for better clarity
metadata <- metadata %>%
  rename(
    'characteristics: patinet ID (Pre=baseline; Post= on treatment)'='patient ID',
    'characteristics: response'= 'response',
    'characteristics: therapy' = 'therapy'
  )

# Create a Seurat object to store the single-cell data and associated metadata
srat <- CreateSeuratObject(counts = counts, meta.data = metadata)

# Display the head and tail of the Seurat object for inspection
head(srat)
tail(srat)


# Display the metadata of the Seurat object
View(srat@meta.data)

# Create a subset of the Seurat object based on the 'therapy' column
filtered_seurat_object <- subset(srat, therapy == "anti-CTLA4+PD1")

# Display the metadata of the filtered Seurat object
View(filtered_seurat_object@meta.data)

# Initialize an empty adjacency matrix
adj.matrix <- NULL

# Display the structure of the filtered Seurat object
str(filtered_seurat_object)

# Extract metadata from the filtered Seurat object
meta <- filtered_seurat_object@meta.data

# Display the dimensions, head, and summary statistics of the metadata
dim(meta)
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

# Calculate the percentage of transcripts mapping to mitochondrial genes
filtered_seurat_object[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

# Display the calculated mitochondrial percentage
filtered_seurat_object[["percent.mt"]]

# Calculate the percentage of transcripts mapping to ribosomal protein genes
filtered_seurat_object[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

# Display the calculated ribosomal protein percentage
filtered_seurat_object[["percent.rb"]]

# Doublet operation (code for doublet operation is omitted)

# Create a violin plot for specified features grouped by 'response'
VlnPlot(filtered_seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1, group.by = "response") +
  theme(plot.title = element_text(size = 10))

# Create scatter plots for various feature combinations
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

# Set up quality control (QC) criteria and create a QC column in metadata
srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass', 'Low_nFeature', srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', srat@meta.data$QC, sep = ','), srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass', 'High_MT', srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT', paste('High_MT', srat@meta.data$QC, sep = ','), srat@meta.data$QC)

# Display the table of QC values
table(srat[['QC']])



# Create logical columns (QC1, QC2, QC3, QC4) based on specified conditions
filtered_seurat_object[['QC1']] <- ifelse((filtered_seurat_object@meta.data$nFeature_RNA > 7500 & filtered_seurat_object@meta.data$response == "Non-responder") | (filtered_seurat_object@meta.data$nFeature_RNA > 4000 & filtered_seurat_object@meta.data$response == "Responder"),"true","false")
filtered_seurat_object[['QC2']] <- ifelse(((filtered_seurat_object@meta.data$nCount_RNA > 38000 | filtered_seurat_object@meta.data$nCount_RNA < 8000) & filtered_seurat_object@meta.data$response == "Non-responder") | (filtered_seurat_object@meta.data$nCount_RNA > 25000 & filtered_seurat_object@meta.data$response == "Responder"),"true","false")
filtered_seurat_object[['QC3']] <- ifelse((filtered_seurat_object@meta.data$percent.mt > 1.67 & filtered_seurat_object@meta.data$response == "Non-responder") | ((filtered_seurat_object@meta.data$percent.mt > 2 | filtered_seurat_object@meta.data$percent.mt <0.37) & filtered_seurat_object@meta.data$response == "Responder"),"true","false")
filtered_seurat_object[['QC4']] <- ifelse(((filtered_seurat_object@meta.data$percent.rb > 10 | filtered_seurat_object@meta.data$percent.rb < 2.31) & filtered_seurat_object@meta.data$response == "Non-responder") | ((filtered_seurat_object@meta.data$percent.rb > 12 | filtered_seurat_object@meta.data$percent.rb <2.85) & filtered_seurat_object@meta.data$response == "Responder"),"true","false")

# Convert logical columns to boolean values
filtered_seurat_object$QC1 <- as.logical(filtered_seurat_object$QC1)
filtered_seurat_object$QC2 <- as.logical(filtered_seurat_object$QC2)
filtered_seurat_object$QC3 <- as.logical(filtered_seurat_object$QC3)
filtered_seurat_object$QC4 <- as.logical(filtered_seurat_object$QC4)

# Combine QC criteria using logical OR
filtered_seurat_object$QC <- filtered_seurat_object$QC1 | filtered_seurat_object$QC2 | filtered_seurat_object$QC3 | filtered_seurat_object$QC4

# Create a violin plot excluding cells failing QC
VlnPlot(subset(filtered_seurat_object, subset = QC == 'FALSE'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1, group.by = "response") +
  theme(plot.title = element_text(size=10))

# Create a new Seurat object (seurat_cleaned) with cells failing QC removed
seurat_cleaned <- subset(filtered_seurat_object, subset = QC == 'FALSE')
seurat_cleaned$combined_group <- paste(seurat_cleaned$response, seurat_cleaned$orig.ident, sep = "")

# View metadata of the filtered Seurat object
View(filtered_seurat_object@meta.data)

# Normalize data and find variable features
seurat_cleaned <- NormalizeData(seurat_cleaned)
seurat_cleaned <- FindVariableFeatures(seurat_cleaned, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_cleaned), 10)

# Plot variable features with labeled top 10 features
plot1 <- VariableFeaturePlot(seurat_cleaned)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Scale data, run PCA, and visualize dimension loadings
all.genes <- rownames(seurat_cleaned)
seurat_cleaned <- ScaleData(seurat_cleaned, features = all.genes)
seurat_cleaned <- RunPCA(seurat_cleaned, features = VariableFeatures(object = seurat_cleaned))
VizDimLoadings(seurat_cleaned, dims = 1:9, reduction = "pca") + 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))

# Create dimensionality heatmap
DimHeatmap(seurat_cleaned, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)

# Plot PCA dimensions
DimPlot(seurat_cleaned, reduction = "pca")

# Plot an elbow plot
ElbowPlot(seurat_cleaned)

# Find neighbors, clusters, and run UMAP
seurat_cleaned <- FindNeighbors(seurat_cleaned, dims = 1:10)
seurat_cleaned <- FindClusters(seurat_cleaned, resolution = 0.5)
seurat_cleaned <- RunUMAP(seurat_cleaned, dims = 1:10, verbose = F)

# Display the table of Seurat clusters
table(seurat_cleaned@meta.data$seurat_clusters)

# Plot UMAP dimensions with labels
DimPlot(seurat_cleaned, label.size = 4, repel = TRUE, label = TRUE)

# Install required packages
install.packages('BiocManager')
install.packages('metap')
BiocManager::install('multtest')

# Extract Seurat clusters for further analysis
seurat_clusters = seurat_cleaned@meta.data$seurat_clusters


# Find conserved markers for cluster 0
cluster0_conserved_markers <- FindConservedMarkers(seurat_cleaned,
                                                   ident.1 = 0,
                                                   grouping.var = "response",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
print(cluster0_conserved_markers)

# Filter conserved markers for cluster 0 based on percentage difference
marqueurs_filtrés <- cluster0_conserved_markers[
  (cluster0_conserved_markers[["Non-responder_pct.1"]] - 
     cluster0_conserved_markers[["Non-responder_pct.2"]]) > 0.8
]
print(marqueurs_filtrés)

# Repeat the process for clusters 1, 2, and 3
cluster1_conserved_markers <- FindConservedMarkers(seurat_cleaned,
                                                   ident.1 = 1,
                                                   grouping.var = "response",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
print(cluster1_conserved_markers)

cluster2_conserved_markers <- FindConservedMarkers(seurat_cleaned,
                                                   ident.1 = 2,
                                                   grouping.var = "response",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
print(cluster2_conserved_markers)

cluster3_conserved_markers <- FindConservedMarkers(seurat_cleaned,
                                                   ident.1 = 3,
                                                   grouping.var = "response",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
print(cluster3_conserved_markers)

# Find conserved markers for cluster 0 based on a different grouping variable (combined_group)
all.conserved_markers_c0 <- FindConservedMarkers(seurat_cleaned,
                                                 ident.1 = 0,
                                                 grouping.var = "combined_group",
                                                 only.pos = TRUE,
                                                 min.diff.pct = 0.25,
                                                 min.pct = 0.25,
                                                 logfc.threshold = 0.25)

# Display information about conserved markers for cluster 0
dim(all.conserved_markers_c0)
table(all.conserved_markers_c0$cluster)
print(all.conserved_markers_c0)

# Select the top 3 markers for each group in the "Non-responder_Pre_pct.1" column
top3_markers <- as.data.frame(all.conserved_markers_c0 %>% group_by("Non-responder_Pre_pct.1") %>% top_n(n = 3, wt = "Non-responder_Pre_pct.1"))


# Find conserved markers for cluster 1, 2, and 3 based on the "combined_group" variable
FindConservedMarkers(seurat_cleaned,
                     ident.1 = 1,
                     grouping.var = "combined_group",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

FindConservedMarkers(seurat_cleaned,
                     ident.1 = 2,
                     grouping.var = "combined_group",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

FindConservedMarkers(seurat_cleaned,
                     ident.1 = 3,
                     grouping.var = "combined_group",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

# Plot features "percent.mt" and "nFeature_RNA"
FeaturePlot(seurat_cleaned, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(seurat_cleaned, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

# Use cell cycle genes and score the cell cycle phases
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seurat_cleaned <- CellCycleScoring(seurat_cleaned, s.features = s.genes, g2m.features = g2m.genes)
table(seurat_cleaned[[]]$Phase)

# Plot features related to mitochondrial genes and ribosomal proteins
FeaturePlot(seurat_cleaned, features = "percent.mt", label.size = 4, repel = T, label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(seurat_cleaned, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(seurat_cleaned, features = "percent.rb", label.size = 4, repel = T, label = T) & theme(plot.title = element_text(size=10))
VlnPlot(seurat_cleaned, features = "percent.rb") & theme(plot.title = element_text(size=10))

# Plot violin plots for RNA count-related features
VlnPlot(seurat_cleaned, features = c("nCount_RNA", "nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))

# Plot features related to cell cycle scoring
FeaturePlot(seurat_cleaned, features = c("S.Score", "G2M.Score"), label.size = 4, repel = T, label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(seurat_cleaned, features = c("S.Score", "G2M.Score")) & 
  theme(plot.title = element_text(size=10))




library(htmltools)
library(glmGamPoi)

# Apply SCTransform with glmGamPoi method, specifying variables to regress
seurat_cleaned <- SCTransform(seurat_cleaned, method = "glmGamPoi", ncells = 8824, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = F)
seurat_cleaned

# Run PCA, UMAP, Find Neighbors, and Find Clusters
seurat_cleaned <- RunPCA(seurat_cleaned, verbose = F)
seurat_cleaned <- RunUMAP(seurat_cleaned, dims = 1:30, verbose = F)
seurat_cleaned <- FindNeighbors(seurat_cleaned, dims = 1:30, verbose = F)
seurat_cleaned <- FindClusters(seurat_cleaned, verbose = F, resolution = 0.5)
table(seurat_cleaned[[]]$seurat_clusters)

# Plot dimensionality reduction
DimPlot(seurat_cleaned, label = T)

# FeaturePlots for specific genes
FeaturePlot(seurat_cleaned, "HLA-DRA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("HLA-DRA")

FeaturePlot(seurat_cleaned, "RGS1") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("RGS1")

FeaturePlot(seurat_cleaned, "IL32") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("IL32")

FeaturePlot(seurat_cleaned, "TRBC2") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("TRBC2")

# Set default assay to "RNA" and normalize data
DefaultAssay(seurat_cleaned) <- "RNA"
seurat_cleaned <- NormalizeData(seurat_cleaned)
seurat_cleaned <- FindVariableFeatures(seurat_cleaned, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_cleaned)
seurat_cleaned <- ScaleData(seurat_cleaned, features = all.genes)

# Find all markers
all.markers <- FindAllMarkers(seurat_cleaned, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# Display dimensionality of all markers
dim(all.markers)

# Display table of markers per cluster
table(all.markers$cluster)

# Select top 3 markers per cluster
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

# Load MonacoImmuneData from celldex package
monaco.ref <- celldex::MonacoImmuneData()

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(DietSeurat(seurat_cleaned))
sce

# SingleR analysis
monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)

# Display tables of pruned labels
table(monaco.main$pruned.labels)
table(monaco.fine$pruned.labels)

# Set identities and visualize
Idents(seurat_cleaned) = "RNA_snn_res.0.5"
seurat_cleaned@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_cleaned@meta.data$monaco.fine <- monaco.fine$pruned.labels
seurat_cleaned <- SetIdent(seurat_cleaned, value = "monaco.main")
DimPlot(seurat_cleaned, label = T, repel = T, label.size = 3) + NoLegend()
DimPlot(seurat_cleaned, group.by = "response") 
DimPlot(seurat_cleaned, group.by = "orig.ident") 

# FeaturePlots for specific genes
FeaturePlot(seurat_cleaned, "HLA-DRA") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(seurat_cleaned, "RGS1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(seurat_cleaned, "IL32") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(seurat_cleaned, "TRBC2") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# Save the Seurat object
saveRDS(seurat_cleaned, file = "srat_tutorial.rds")

# Rename identities
seurat_cleaned <- RenameIdents(seurat_cleaned, 'Monocytes' = '0')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'CD4+ T cells' = '1')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'B cells' = '2')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'T cells' = '3')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'CD8+ T cells' = '4')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'DendriticI cells' = '5')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'NK cells' = '6')
seurat_cleaned <- RenameIdents(seurat_cleaned, 'Neutrophils' = '7')

# Find markers for cluster 0
cluster0.markers <- FindMarkers(seurat_cleaned, ident.1 = 0)
head(cluster0.markers)

# Find all markers and filter based on log2FC
seurat_cleaned.markers <- FindAllMarkers(seurat_cleaned, only.pos = TRUE)
seurat_cleaned.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Find markers using ROC test for cluster 0
cluster0.markers <- FindMarkers(seurat_cleaned, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(seurat_cleaned, features = c("MS4A1", "CD79A"), slot = "counts")

# FeaturePlot for specific genes
FeaturePlot(seurat_cleaned, features = c("CD14", "SERPINA1", "LYZ", "RP11-1143G9.4", "LRP1", "TYROBP", "CST3", "CLEC7A", "AIF1"))

# Select top 10 variable features and create a heatmap
top10 <- head(VariableFeatures(seurat_cleaned), 10)

seurat_cleaned.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(seurat_cleaned, features = top10$gene) + NoLegend()

# Rename identities for better interpretation
seurat_cleaned <- RenameIdents(seurat_cleaned, '0' = 'Macrophages')
seurat_cleaned <- RenameIdents(seurat_cleaned, '1' = 'T cells')
seurat_cleaned <- RenameIdents(seurat_cleaned, '2' = 'B cells')
seurat_cleaned <- RenameIdents(seurat_cleaned, '3' = 'Gamma delta T cells')
seurat_cleaned <- RenameIdents(seurat_cleaned, '5' = 'Dendritic cells')
seurat_cleaned <- RenameIdents(seurat_cleaned, '4' = 'NK cells2')
seurat_cleaned <- RenameIdents(seurat_cleaned, '6' = 'NK cells1')
seurat_cleaned <- RenameIdents(seurat_cleaned, '7' = 'Neutrophils')

# Plotting the UMAP with labeled identities
DimPlot(seurat_cleaned, label = T , repel = T, label.size = 3) + NoLegend()

# Find conserved markers for B cells
FindConservedMarkers(seurat_cleaned,
                     ident.1 = "B cells",
                     grouping.var = "response",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

# VlnPlot for MS4A1 expression
VlnPlot(seurat_cleaned, features = c("MS4A1"), slot = "counts")

# UMAP with labeled combined groups
DimPlot(seurat_cleaned, reduction = "umap", split.by = "combined_group")

# Find conserved markers for B cells and plot them
bcells.markers <- FindConservedMarkers(seurat_cleaned, ident.1 = "B cells", grouping.var = "response", verbose = FALSE)
head(bcells.markers)

# FeaturePlot for MS4A1, CD19, and LYZ
FeaturePlot(seurat_cleaned, features = c("MS4A1", "CD19", "LYZ"), split.by = "combined_group", max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")

# VlnPlots for specific markers
plots <- VlnPlot(seurat_cleaned, features = c("MS4A1", "CD19", "LYZ"), split.by = "combined_group", group.by = "monaco.main",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# VlnPlots for CD14 and MS4A1
plots <- VlnPlot(seurat_cleaned, features = c("CD14", "MS4A1"), split.by = "combined_group", group.by = "monaco.main",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# FeaturePlot for CD14 and MS4A1
FeaturePlot(seurat_cleaned, features = c("CD14", "MS4A1"), split.by = "combined_group", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")

# VlnPlots for CCL5 and CD4
plots <- VlnPlot(seurat_cleaned, features = c("CCL5", "CD4"), split.by = "combined_group", group.by = "monaco.main",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# FeaturePlot for CCL5 and CD4
FeaturePlot(seurat_cleaned, features = c("CCL5", "CD4"), split.by = "combined_group", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")

# VlnPlots for CD28
plots <- VlnPlot(seurat_cleaned, features = c("CD28"), split.by = "combined_group", group.by = "monaco.main",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# FeaturePlot for CD28
FeaturePlot(seurat_cleaned, features = c("CD28"), split.by = "combined_group", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")

# Extract meta data from seurat_cleaned
meta_data <- seurat_cleaned@meta.data

# Count occurrences of combinations in seurat_cleaned
table_count_seurat <- table(meta_data$orig.ident, meta_data$response)

# Assuming 'filtered_seurat_object' is another Seurat object
sample <- filtered_seurat_object

# Perform necessary pre-processing steps on 'sample'
sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = 0.5)
sample <- RunUMAP(sample, dims = 1:10, verbose = F)

# Extract meta data from the processed 'sample'
meta_data_sample <- sample@meta.data

# Count occurrences of combinations in the processed 'sample'
table_count_sample <- table(meta_data_sample$orig.ident, meta_data_sample$response)

# Create a DimPlot to visualize the UMAP with color-coded clusters
DimPlot(sample, reduction = "umap", split.by = "combined_group")
