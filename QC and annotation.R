library(SoupX)
library(ggsci)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(harmony)

setwd("")

######SOUPX load data (HB01 as one example)
HB01 = load10X("HB01/")
HB01 = autoEstCont(HB01)
HB01 = adjustCounts(HB01)

HB01 <- CreateSeuratObject(counts = HB01, project = 'HB01', min.cells = 3,min.features = 200)
ncol(HB01)
HB01[["percent.mt"]] <- PercentageFeatureSet(HB01, pattern = "^MT-")
HB01<- subset(HB01, subset = nFeature_RNA > 200 & percent.mt < 25)
ncol(HB01)
HB01 <- HB01[!grepl("^MT-", rownames(HB01)), ]
HB01 <- HB01[!grepl('^RP[SL]', rownames(HB01)), ]
HB01 <- NormalizeData(HB01, normalization.method = "LogNormalize", scale.factor = 10000)
HB01 <- FindVariableFeatures(HB01, selection.method = "vst", nfeatures = 2000)
HB01 <- ScaleData(HB01)
HB01 <- RunPCA(HB01, features = VariableFeatures(object = HB01),nps=20,verbose = TRUE)
HB01 <- FindNeighbors(HB01,reduction="pca",dims = 1:15)
HB01 <- FindClusters(HB01, resolution = 0.1)
HB01<-RunUMAP(HB01,dims = 1:20,reduction = "pca")

DimPlot(HB01, reduction = "umap",seed = 3,label = TRUE, label.size = 5)+scale_color_igv("default")

####Find doublets####
sweep.res.list_HB01 <- paramSweep_v3(HB01, PCs = 1:40, sct = FALSE)
sweep.stats_HB01 <- summarizeSweep(sweep.res.list_HB01, GT = FALSE)
bcmvn_HB01 <- find.pK(sweep.stats_HB01)
opt_pK <- as.numeric(as.vector(bcmvn_HB01$pK[which.max(bcmvn_HB01$BCmetric)]))
print(opt_pK)

annotations <- HB01@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.075*ncol(HB01@assays$RNA@data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HB01 <- doubletFinder_v3(HB01, 
                         PCs = 1:10, 
                         pN = 0.25, 
                         pK = opt_pK, 
                         nExp = nExp_poi, 
                         reuse.pANN = FALSE, 
                         sct = FALSE)
colnames(HB01@meta.data)[ncol(HB01@meta.data)]="DoubletFinder"
table(HB01$DoubletFinder)
DimPlot(HB01,reduction = "umap",pt.size = 0.1,group.by = "DoubletFinder",seed = 3)

####merge all the data as 'HB'

HB<- NormalizeData(HB, normalization.method = "LogNormalize", scale.factor = 10000)
HB<- FindVariableFeatures(HB, selection.method = "vst", nfeatures = 1000)
HB<- ScaleData(HB, features = VariableFeatures(HB))
HB <- RunPCA(HB, features = VariableFeatures(object = HB),nps=50,verbose = TRUE)
ElbowPlot(HB,ndims = 50)
HB <- FindNeighbors(HB,reduction = "pca",dims = 1:25)
HB<-RunUMAP(HB,dims = 1:25,reduction = "pca")
HB <- FindClusters(HB, resolution = 0.2)
DimPlot(HB, reduction = "umap",label = TRUE, label.size = 5,raster = FALSE)+scale_color_d3("category20c")
DimPlot(HB, reduction = "umap",seed = 3,group.by = "orig.ident")+scale_color_igv("default")

HB<- subset(HB,subset = DoubletFinder =="Singlet")

HB@meta.data$sample_type[HB@meta.data$orig.ident%in%c("HB01","HB03","HB04","HB05","HB06T",
                                                      "HB07T","HB09T","HB10T","HB11","HB08",
                                                      "HB14T","HB16T","HB17T")]<-"tumor"

HB@meta.data$sample_type[HB@meta.data$orig.ident%in%c("HB06P","HB07P","HB09P","HB10P","HB14P",
                                                      "HB16P","HB17P")]<-"non_tumor"
table(HB$sample_type)
DimPlot(HB, reduction = "umap",group.by = "sample_type",raster=FALSE,label = T)+
  scale_color_d3("category20c")

###First-step annotation

FeaturePlot(HB,features=c("ALB","IGF2","H19","PTPRC","CD3E","NKG7","CD4","CD8A","GZMK","IFNG","FOXP3","NCAM1","CD1C",
                          "CLEC9A","CD68","APOE","TPSAB1","LILRA4","CD79A",
                          "AHSP","GYPA","ACTA2","PDGFRA","MCAM","ELANE","CD34",
                          "SPINK2","JCHAIN","PLVAP","CLEC4M","MKI67","EPCAM"),ncol = 6,raster=FALSE)

HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("1","2","5","15")]<-"Malignancy"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("12")]<-"Hepatocyte"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("0","11")]<-"NK/T"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("13")]<-"Cycling"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("4","7","10","16")]<-"Myeloid"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("3")]<-"B_cell"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("9")]<-"EPC"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("6")]<-"Endothelial_cell"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("8")]<-"Fibroblast"
HB@meta.data$celltype[HB@meta.data$seurat_clusters%in%c("14")]<-"Cholangiocyte"

DimPlot(HB, reduction = "umap",group.by = "celltype",label =T,label.size = 5)+
  scale_color_d3("category20c")+theme(legend.position = "bottom")


####Sub-clustering (NK/T cell as one example)
nkt <- subset(HB,idents="NK/T")

nkt<- FindVariableFeatures(nkt, selection.method = "vst", nfeatures = 500)
nkt<- ScaleData(nkt, features = VariableFeatures(nkt))
nkt <- RunPCA(nkt, features = VariableFeatures(object = nkt),nps=50,verbose = TRUE)
nkt <- nkt %>% RunHarmony("orig.ident",plot_convergence=TRUE)
nkt<-RunUMAP(nkt,dims = 1:30,reduction = "harmony")
nkt <- FindNeighbors(nkt,reduction = "harmony",dims = 1:30)
nkt <- FindClusters(nkt, resolution = 0.4)

# Visualization
DimPlot(nkt, reduction = "umap", label = TRUE)+scale_color_d3("category20c")
DimPlot(nkt, reduction = "umap", group.by = "orig.ident", repel = TRUE)+scale_color_igv("default")+NoLegend()
DimPlot(nkt, reduction = "umap", group.by = "sample_type", repel = TRUE)+scale_color_d3("category20c")
#
nkt$annotation <- "T_cell"
nkt$annotation[nkt$seurat_clusters%in%c('1')] <- "NK_CX3CR1"
nkt$annotation[nkt$seurat_clusters%in%c('0','10')] <- "NK_NCAM1"
nkt$annotation[nkt$seurat_clusters%in%c('5')] <- "gdT"
nkt$annotation[nkt$seurat_clusters%in%c('4')] <- "MAIT"

##
Idents(nkt) <- "annotation"
t_cell <- subset(nkt,idents = c("T_cell"))

t_cell<- FindVariableFeatures(t_cell, selection.method = "vst", nfeatures = 500)
t_cell<- ScaleData(t_cell, features = VariableFeatures(t_cell))
t_cell <- RunPCA(t_cell, features = VariableFeatures(object = t_cell),nps=50,verbose = TRUE)
t_cell <- t_cell %>% RunHarmony("orig.ident",plot_convergence=TRUE)
t_cell<-RunUMAP(t_cell,dims = 1:30,reduction = "harmony")
t_cell <- FindNeighbors(t_cell,reduction = "harmony",dims = 1:30)
t_cell <- FindClusters(t_cell, resolution = 1)
DimPlot(t_cell, reduction = "umap", label = TRUE)+scale_color_igv("default")

t_cell$annotation <- "na"
t_cell$annotation[t_cell$seurat_clusters%in%c("0","11")] <- "CD8_c1_GZMK"
t_cell$annotation[t_cell$seurat_clusters=="1"] <- "CD4_c1_CD40LG"
t_cell$annotation[t_cell$seurat_clusters%in%c("2","16")] <- "CD8_c2_ZNF683"
t_cell$annotation[t_cell$seurat_clusters%in%c("3")] <- "Naive_T"
t_cell$annotation[t_cell$seurat_clusters%in%c("4","10")] <- "CD8_c3_CX3CR1"
t_cell$annotation[t_cell$seurat_clusters=="5"] <- "CD4_c2_FOXP3"
t_cell$annotation[t_cell$seurat_clusters%in%c("7","20")] <- "CD8_c4_CXCL13"
t_cell$annotation[t_cell$seurat_clusters=="8"] <- "NK_c3_TYROBP"
t_cell$annotation[t_cell$seurat_clusters=="9"] <- "CD4_c3_IL7R"
t_cell$annotation[t_cell$seurat_clusters=="12"] <- "Cycling_T"
t_cell$annotation[t_cell$seurat_clusters=="13"] <- "ILC3_RORC"
t_cell$annotation[t_cell$seurat_clusters=="17"] <- "ILC1_TBX21"
t_cell$annotation[t_cell$seurat_clusters=="18"] <- "ILC2_GATA3"

DimPlot(t_cell, reduction = "umap",group.by = "annotation", label = TRUE)+scale_color_igv("default")