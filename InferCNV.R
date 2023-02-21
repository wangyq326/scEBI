library(rjags)
library(infercnv)
library(ggsci)
library(Seurat)
library(ggplot2)
library(dplyr)
library(usethis)
library(devtools)
library(AnnoProbe)

###TEST#######
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)


#####
setwd(" ")
HB_tumor <- readRDS("tumor annotation.rds")
Idents(HB_tumor) <- "orig.ident"

immune=readRDS("immune annotation.rds")
Idents(immune)="annotation"
table(immune$annotation)
immune=subset(immune,idents=c("Macrophage","DC","T_cell"))
infer_one=merge(HB_tumor,immune)
infer_one$group=infer_one$orig.ident
infer_one$group[infer_one$annotation%in%c("Macrophage","DC","T_cell")]="ref"
table(infer_one$group)
Idents(infer_one) <- "group"
###infercnv
library(Seurat)
library(magrittr)
library(ggplot2)
infer_one_matrix<-GetAssayData(infer_one@assays$RNA, slot="counts")
infer_one_annotation<-data.frame(rownames(infer_one@meta.data),infer_one@meta.data$group)
write.table(infer_one_annotation,"/infercnv/infer_one_annotation.txt",sep = '\t', quote = F,col.names = F,row.names = F)

infer_one_obj = CreateInfercnvObject(raw_counts_matrix=infer_one_matrix,
                                     annotations_file="/infercnv/infer_one_annotation.txt",
                                     delim="\t",
                                     gene_order_file="reference/gene-locate-board.txt",
                                     ref_group_names=c("ref")) 

infer_one_obj = infercnv::run(infer_one_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir="/infercnv",
                              cluster_by_groups=TRUE, 
                              denoise=TRUE,   
                              HMM=FALSE)  