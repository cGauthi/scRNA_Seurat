BiocManager::install('limma')

library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(openxlsx)
library(tidyverse)
library(limma)


################################################  CLUSTER MARKER VISUALIZATION ###############################

# Load Seurat data set
scrna_filtered <- readRDS("mono_filtered_SeuratObject.rds")

# Create cell cluster UMAP
DimPlot(scrna_filtered, reduction = 'umap', label = TRUE) + NoLegend()

# Create HTO-based UMAPs
DimPlot(scrna_filtered, group.by = 'HTO_classification')
DimPlot(scrna_filtered, group.by = 'HTO_classification', split.by = 'HTO_classification', ncol = 4)

# Create marker lists of interest
classic_marks <- c('CCR2', 'CD14', 'CD163', 'CD36', 'SELL', 'SERPINB2', 'CLEC4E', 'HLA-DRB1', 'SIRPA', 'TREM1')
int_marks <- c('CD40', 'CD74', 'HLA-DRA', 'CCR5', 'ITGAE', 'ITGAX', 'CSF1R', 'IRF8', 'RELA', 'TLR7', 'MARCO')
nonclass_marks <- c('FCGR3A', 'ITGAL', 'SIGLEC10', 'MSR1', 'ADAM17', 'C1QA', 'LTB', 'MICA', 'TLR9', 'TNF')
common_marks <- c('TLR2', 'ITGB2', 'ITGAM', 'CTSD', 'CTSA', 'NLRP3', 'CLEC7A', 'BST1', 'STAB1', 'IRAK3')

## MARKER IDENTIFICATION ##
# Listing cluster markers
# c1_markers <- FindMarkers(scrna_filtered, ident.1 = 1, min.pct = 0.25)
# head(c1_markers, n = 10)

# # Listing differences betwen individual clusters
DEGenes_xGroup <- FindMarkers(scrna_filtered, group.by = 'HTO_classification', ident.1 = c('AD1', 'AD2', 'AD3', 'AD4'), ident.2 = c('Control1', 'Control2'))
# head(DEmarks, n = 10)
write.xlsx(DEGenes_xGroup, 'DEGenes_xGroup_reducedControls.xlsx', rowNames = TRUE)

# List all positive markers for each cluster
scrna_filtered_markers <- FindAllMarkers(scrna_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scrna_filtered_markers <- scrna_filtered_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
write.xlsx(scrna_filtered_markers, 'top_markers_dend.xlsx', rowNames = TRUE)

## MARKER VISUALIZATION ##
# Violin plots for expression levels by cluster
VlnPlot(scrna_filtered, features = classic_marks, pt.size = 0.2, ncol = 5)
VlnPlot(scrna_filtered, features = rownames(DEGenes_xGroup), pt.size = 0.2, ncol = 5, group.by = 'HTO_classification')

# Feature map of marker genes on the sctransform embedding.
FeaturePlot(scrna_filtered, features =  rownames(DEGenes_xGroup), pt.size = 0.2, ncol = 4)

# Heatmap showing top genes for each cluster (requires a full model - sct is a reduced model)
#DoHeatmap(scrna_filtered, features = scrna_filtered_markers$gene) + NoLegend()

# Assignment of cll type to clusters via canonical markers (fill new_cluster_ids with celltypes in order)
new_cluster_ids <- c('Mono1', 'Mono2', 'Mono2', 'Interferon-induced', 'Mono2', 'CD14+/CD16+', 'Mono2', 'Mono1', 'CD14-/CD16+', 'Mono2', 'Mono1', 'Cytotoxic')
names(new_cluster_ids) <- levels(scrna_filtered)
scrna_filtered <- RenameIdents(scrna_filtered, new_cluster_ids)
DimPlot(scrna_filtered, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

