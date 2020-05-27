# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# https://satijalab.org/seurat/v3.0/sctransform_vignette.html
# https://satijalab.org/seurat/v3.1/integration.html

library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
options(future.globals.maxSize = 4000 * 1024^2)

##################################  DATA PREPARATION  ############################################
# Read in the data
data <- Read10X(data.dir = '../../sequencing_data/scRNA/684/filtered_feature_bc_matrix')


# Extract UMI and Hashtag Oligo matrices
data_UMI <- data$`Gene Expression`
data_HTO <- as.matrix(data$`Antibody Capture`)

joint_barcodes <- intersect(colnames(data_UMI), colnames(data_HTO))

# Subset UMI and Hashtag Oligo matrices by shared barcodes
data_UMI <- data_UMI[, joint_barcodes]
data_HTO <- as.matrix(data_HTO[, joint_barcodes])

rownames(data_HTO)


######################### MITOCHONDRIAL REGRESSION ############################################
# Create Seurat object
scrna_base <- CreateSeuratObject(counts = data_UMI, project = 'BRI-684')

# Store mitochondrial gene statatistics in your Seurat object
scrna_base[['percent_mt']] <- PercentageFeatureSet(scrna_base, pattern = '^MT-')

# Check mitochondrial metrics (can use to set filter limits at key points)
VlnPlot(scrna_base, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
scat1 <- FeatureScatter(scrna_base, feature1 = 'nCount_RNA', feature2 = 'percent_mt')
scat2 <- FeatureScatter(scrna_base, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scat1 + scat2

# # SCTransform - new command that combines NormalizeData, ScaleData, and FindVariableFeatures - MUCH FASTER (and better), but loses some ability to process non-DE genes (i.e. heatmapping)
scrna_base <- SCTransform(scrna_base, vars.to.regress = 'percent_mt', verbose = TRUE)

############################  DEMULTIPLEXING ###########################################
# Add Hashtag Oligo data as an independent assay to test with
scrna_base[["HTO"]] <- CreateAssayObject(counts = data_HTO)

# Normalize data using HTOs and centered log-ratio (CLR) transformation
scrna_base <- NormalizeData(scrna_base, assay = "HTO", normalization.method = 'CLR')

# Demultiplex cells based on their HTO levels
scrna_base <- HTODemux(scrna_base, assay = 'HTO', positive.quantile = 0.99)

# Global classification results (Singlets/Doublets/Negatives)
table(scrna_base$HTO_classification.global)
Idents(scrna_base) <- 'HTO_classification.global'
VlnPlot(scrna_base, features = 'nCount_RNA', pt.size = 0.1, log = TRUE)

# Enrichment plots by HTO group
Idents(scrna_base) <- 'HTO_maxID'
RidgePlot(scrna_base, assay = 'HTO', features = rownames(scrna_base[['HTO']])[1:8], ncol = 2)

# Compare HTO signals between groups for confirmation using scatterplot (optional)
Idents(scrna_base) <- 'HTO_maxID'
FeatureScatter(scrna_base, feature1 = 'hto_Hashtag3', feature2 = 'hto_Hashtag7')

# tSNE Plots of HTO groups (will take a few minutes)
Idents(scrna_base) <- 'HTO_classification.global'
scrna_subset <- subset(scrna_base, idents = 'Negative', invert = TRUE)                     # Removal of Negative cell
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = scrna_subset, assay = 'HTO'))))   # Calculate a distance matrix for HTOs
scrna_subset <- RunTSNE(scrna_subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(scrna_subset)

# HTO-based Heatmap (shows singlet and doublet grouping)
HTOHeatmap(scrna_subset, assay = 'HTO', ncells = 5000)                     # ncells is the number of cells to subset


###########################   BATCH EFFECT REMOVAL #######################
# Select only Singlet cells and split into difference batches based on their HTO
scrna_split <- subset(scrna_base, idents = 'Singlet')
scrna_split <- SplitObject(scrna_split, split.by = 'HTO_classification')

# Run SCTransform on all batches
for (i in 1:length(scrna_split)) {
  scrna_split[[i]] <- SCTransform(scrna_split[[i]], verbose = TRUE)
}

# Determine and prepare top shared features to be used as anchors for data integration
scrna_features <- SelectIntegrationFeatures(object.list = scrna_split, nfeatures = 3000)
scrna_split <- PrepSCTIntegration(object.list = scrna_split, anchor.features = scrna_features, verbose = TRUE)

# Set anchors between batches and reintegrate into one Seurat dataset
scrna_anchors <- FindIntegrationAnchors(object.list = scrna_split, normalization.method = 'SCT', anchor.features = scrna_features, verbose = TRUE)
scrna_integrated <- IntegrateData(anchorset = scrna_anchors, normalization.method = 'SCT', verbose = TRUE)

# Determine principal components and Visualize debatched and integrated dataset
scrna_integrated <- RunPCA(scrna_integrated, verbose = FALSE)
scrna_integrated <- RunUMAP(scrna_integrated, dims = 1:30)

DimPlot(scrna_integrated, group.by = 'HTO_classification')

################################# CELLTYPE IDENTIFICATION  #############################################3
# Further filtering of unwanted cells (low feature count and high mitochondrial percentage) - may have LARGE IMPACT - Use good judgment
scrna_filtered <- subset(scrna_integrated, subset = nFeature_RNA > 200 & percent_mt < 20)

# Select out the top 1000 most variable features/genes
scrna_filtered <- SCTransform(scrna_filtered)

# Run PCA on filtered dataset
scrna_filtered <- RunPCA(scrna_filtered, features = VariableFeatures(scrna_filtered))

# Determine the dimensionality of the data through use of JackStraw plots (advanced) or Elbow plots (quick and dirty)
# JackStraw (computation heavy) can only be run on non-SCT data.  Use normal pipeline (Normalize, FindVariable, Scale) to use this rather than Elbow
#scrna_singlet <- JackStraw(scrna_singlet, num.replicate = 100)
#scrna_singlet <- ScoreJackStraw(scrna_singlet, dims = 1:20)
#JackStrawPlot(scrna_singlet, dims = 1:20)

# ElbowPlot (quick approximation)
ElbowPlot(scrna_filtered)

# Select appropriate PCs for cluster and UMAP based on PCElbowPlot
scrna_filtered <- FindNeighbors(scrna_filtered, reduction = 'pca', dims = 1:8)
scrna_filtered <- FindClusters(scrna_filtered, resolution = 0.6, verbose = TRUE)
scrna_filtered <- RunUMAP(scrna_filtered, reduction = 'pca', dims = 1:10)

# Project singlet identities on UMAP plot
DimPlot(scrna_filtered, group.by = 'HTO_classification')


########################  DIFFERENTIAL FEATURE PROCESSING ##########################################
top10 <- head(VariableFeatures(scrna_filtered), 10)

# Visualized variable features.  If run on a nonSCT model, there will be more non-variables visualized due to inclusion in the model.
featplot1 <- VariableFeaturePlot(scrna_filtered)
featplot2 <- LabelPoints(plot = featplot1, points = top10, repel = TRUE)
featplot2

# Data scaling (sets mean expression to 0 and scales expression for total variance across cells = 1.  Prevents domination by highly-expressed genes)
#all_genes <- rownames(scrna_filtered)
#scrna_filtered <- ScaleData(scrna_filtered, features = all_genes, vars.to.regress = 'percent_mt')           # can be run only on variable genes by leaving out features (default) - ruins heatmapping

################################# DIMENSIONALITY VISUALIZATION ############################################
scrna_filtered <- RunPCA(scrna_filtered, features = VariableFeatures(object = scrna_filtered))

# Visualize the principal components (multiple options)
print(scrna_filtered[['pca']], dims = 1:10, nfeatures = 5)
VizDimLoadings(scrna_filtered, dims = 1:10, reduction = 'pca')
DimPlot(scrna_filtered, reduction = 'pca')
DimHeatmap(scrna_filtered, dims = 1:10, cells = 500, balanced = TRUE)


#####################################  CELL CLUSTERING VISUALIZATION ###################################################
# For SCTransform Method (still better)
# Construct KNN graph based on dimensions chosen
scrna_filtered <- FindNeighbors(scrna_filtered, dims = 1:10, verbose = FALSE)

# Use Louvain algorith or SLM to iteratively group cells.  Resolution parameter is variable, but 0.4-1.2 works well for datasets w/ around 3000 cells
scrna_filtered <- FindClusters(scrna_filtered, resolution = 0.6, verbose = FALSE)

# Run UMAP (better than tSNE)
scrna_filtered <- RunUMAP(scrna_filtered, dims = 1:10, verbose = TRUE)
DimPlot(scrna_filtered, reduction = 'umap', label = TRUE) + NoLegend()


################################################  CLUSTER MARKER VISUALIZATION ###############################
# Violin plots for expression levels by cluster
VlnPlot(scrna_filtered, features = c('CD19', 'MS4A1', 'CD27', 'CD14', 'LYZ', 'GNLY', 'S100A4', 'FCGR3A'),
      pt.size = 0.2, ncol = 4)

# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(scrna_filtered, features =  c('MS4A1', 'CD14', 'LYZ', 'PPBP', 'CD8A', 'IL7R', 'CCR7', 'S100A4', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP'),
            pt.size = 0.2, ncol = 3)

# Listing cluster markers (hmmm questionable... filter issue?)
c1_markers <- FindMarkers(scrna_filtered, ident.1 = 1, min.pct = 0.25)
head(c1_markers, n = 10)

# Listing differences betwen individual clusters (also questionable)
c5_c0_DEmarks <- FindMarkers(scrna_filtered, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
head(c5_c0_DEmarks, n = 10)

# List all positive markers for each cluster
scrna_filtered_markers <- FindAllMarkers(scrna_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scrna_filtered_markers <- scrna_filtered_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)

# Heatmap showing top genes for each cluster (requires a full model - sct is a reduced model)
#DoHeatmap(scrna_filtered, features = scrna_filtered_markers$gene) + NoLegend()

# Assignment of cll type to clusters via canonical markers (fill new_cluster_ids with celltypes in order)
new_cluster_ids <- c('B1', 'B2', 'B3', 'Memory CD4+ (?)', 'B4', 'Unknown1', 'NK', 'Unknown2', 'CD14+ Monocytes', 'FCGR3A+ Mono')
names(new_cluster_ids) <- levels(scrna_filtered)
scrna_filtered <- RenameIdents(scrna_filtered, new_cluster_ids)
DimPlot(scrna_filtered, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()



