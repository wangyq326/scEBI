##reference:https://cole-trapnell-lab.github.io/monocle3/


library(monocle3)
library(Seurat)
library(ggsci)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(harmony)
rm(list=ls())

nkt <- NormalizeData(nkt, normalization.method = "LogNormalize", scale.factor = 10000)
nkt<- FindVariableFeatures(nkt, selection.method = "vst", nfeatures = 1000)
nkt<- ScaleData(nkt, features = VariableFeatures(nkt))
nkt <- RunPCA(nkt, features = VariableFeatures(object = nkt),nps=50,verbose = TRUE)
nkt <- nkt %>% RunHarmony("orig.ident",plot_convergence=TRUE)
nkt<-RunUMAP(nkt,dims = 1:8,reduction = "harmony")
nkt <- FindNeighbors(nkt,reduction = "harmony",dims = 1:8)
nkt <- FindClusters(nkt, resolution = 0.1)
DimPlot(nkt, reduction = "umap", label = TRUE,group.by = "annotation")+scale_color_d3("category20c")
DimPlot(nkt, reduction = "umap", label = TRUE)+scale_color_d3("category20c")


##
data <- GetAssayData(nkt, assay = 'RNA', slot = 'counts')
cell_metadata <- nkt@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds)
#
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="annotation2") + ggtitle('cds.umap')

##
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(nkt, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="annotation2") + ggtitle('int.umap')

## 
cds <- cluster_cells(cds,reduction_method = "UMAP",k = 13)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = p1+p2
p
##
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds,
           color_cells_by = "annotation",
           label_cell_groups=T,group_label_size = 3,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)+scale_color_d3("category20c")


####progenitors
hb_merge=readRDS("HB_merge_annotation.rds")
Idents(hb_merge) <- "celltype"
hb_merge <- subset(hb_merge,idents=c("Hepatocyte","Endothelial_cell","Fibroblast","Malignancy","Cholangiocyte"),invert=T)
Idents(hb_merge) <- "annotation2"
hb_merge <- subset(hb_merge,idents=c("Cycling immune cell","Plasmablast","Plasma_cell","gdT","T cell","ILC","B cell",
                                     "Pre_B","Pro_B","Pre-pro_B","pDC"),invert=T)

hb_merge <- NormalizeData(hb_merge, normalization.method = "LogNormalize", scale.factor = 10000)
hb_merge<- FindVariableFeatures(hb_merge, selection.method = "vst", nfeatures = 400)
hb_merge<- ScaleData(hb_merge, features = VariableFeatures(hb_merge))
hb_merge <- RunPCA(hb_merge, features = VariableFeatures(object = hb_merge),nps=50,verbose = TRUE)
hb_merge <- hb_merge %>% RunHarmony("orig.ident",plot_convergence=TRUE)
hb_merge<-RunUMAP(hb_merge,dims = 1:8,reduction = "harmony")
hb_merge <- FindNeighbors(hb_merge,reduction = "harmony",dims = 1:8)
hb_merge <- FindClusters(hb_merge, resolution = 0.1)
DimPlot(hb_merge, reduction = "umap", label = TRUE,group.by = "annotation2")+scale_color_d3("category20c")
DimPlot(hb_merge, reduction = "umap", label = TRUE)+scale_color_d3("category20c")

Idents(hb_merge) <- "annotation2"
hb_merge <- subset(hb_merge,idents=c("Mast cell"),invert=T)

##
data <- GetAssayData(hb_merge, assay = 'RNA', slot = 'counts')
cell_metadata <- hb_merge@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds)
#
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="annotation2") + ggtitle('cds.umap')

##
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(hb_merge, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="annotation2") + ggtitle('int.umap')

## 
cds <- cluster_cells(cds,reduction_method = "UMAP",k = 13)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = p1+p2
p
##
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds,
           color_cells_by = "annotation2",
           label_cell_groups=F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)+scale_color_d3("category20c")
