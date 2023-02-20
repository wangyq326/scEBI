suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
###load data
HB01 <- read.csv("TCR/tcr/HB01_filtered_contig_annotations.csv")
HB03 <- read.csv("TCR/tcr/HB03_filtered_contig_annotations.csv")
HB04 <- read.csv("TCR/tcr/HB04_filtered_contig_annotations.csv")
HB05 <- read.csv("TCR/tcr/HB05_filtered_contig_annotations.csv")
HB06T <- read.csv("TCR/tcr/HB06T_filtered_contig_annotations.csv")
HB06P <- read.csv("TCR/tcr/HB06P_filtered_contig_annotations.csv")
HB07T <- read.csv("TCR/tcr/HB07T_filtered_contig_annotations.csv")
HB07P <- read.csv("TCR/tcr/HB07P_filtered_contig_annotations.csv")
HB08 <- read.csv("TCR/tcr/HB08_filtered_contig_annotations.csv")
HB09T <- read.csv("TCR/tcr/HB09T_filtered_contig_annotations.csv")
HB09P <- read.csv("TCR/tcr/HB09P_filtered_contig_annotations.csv")
HB10T <- read.csv("TCR/tcr/HB10T_filtered_contig_annotations.csv")
HB10P <- read.csv("TCR/tcr/HB10P_filtered_contig_annotations.csv")
HB11 <- read.csv("TCR/tcr/HB11T_filtered_contig_annotations.csv")
HB14T <- read.csv("TCR/tcr/HB14T_filtered_contig_annotations.csv")
HB14P <- read.csv("TCR/tcr/HB14P_filtered_contig_annotations.csv")
HB16T <- read.csv("TCR/tcr/HB16T_filtered_contig_annotations.csv")
HB16P <- read.csv("TCR/tcr/HB16P_filtered_contig_annotations.csv")
HB17T <- read.csv("TCR/tcr/HB17T_filtered_contig_annotations.csv")
HB17P <- read.csv("TCR/tcr/HB17P_filtered_contig_annotations.csv")


contig_list <- list(HB01,HB03,HB04,HB05,HB06T,HB06P,HB07T,HB07P,HB08,HB09T,HB09P,HB10T,HB10P,
                    HB11,HB14T,HB14P,HB16T,HB16P,HB17T,HB17P)
combined <- combineTCR(contig_list, 
                       samples = c("HB01", "HB03", "HB04", "HB05", "HB06T","HB06P","HB07T","HB07P","HB08",
                                   "HB09T","HB09P","HB10T","HB10P","HB11","HB14T","HB14P","HB16T","HB16P","HB17T","HB17P"), 
                       ID=c("T","T","T","T","T","P","T","P","T","T","P","T","P","T","T","P","T","P","T","P"),
                       cells ="T-AB", filterMulti = FALSE)

for (i in seq_along(combined)) {
  for (j in 1:nrow(combined[[i]])) {
    combined[[i]][j,1] <- gsub(paste(c(as.character(combined[[i]][j,2]),"_",as.character(combined[[i]][j,2])),collapse = ""),as.character(combined[[i]][j,2]),combined[[i]][j,1])    }
}
head(combined[[1]])

saveRDS(combined,"combined_TCR_nofilter.rds")
combined=readRDS("combined_TCR_nofilter.rds")
#write.csv(cd8@meta.data,"cd8 metadata.csv")##
cd8meta=read.csv("cd8 metadata.csv",header = T)
rownames(cd8meta)=cd8meta$X
cd8meta=cd8meta[,-1]

nkt4=readRDS("nkt_annotation.rds")
Idents(nkt4) <- "annotation"
cd8 <- subset(nkt4,idents=c("CD8_c1_GZMK","CD8_c2_ZNF683","CD8_c3_CX3CR1","CD8_c4_CXCL13"))
cd8@meta.data=cd8meta
CD8_TCR<- combineExpression(combined, cd8, 
                            cloneCall="gene+nt", 
                            cloneTypes=c(No=0,Single=1, Small=5, Medium=10,Large=25,Hyperexpanded=80))

table(CD8_TCR@meta.data$cloneType)

slot(CD8_TCR, "meta.data")$cloneType <- factor(slot(CD8_TCR, "meta.data")$cloneType, 
                                               levels = c("Hyperexpanded (25 < X <= 80)",
                                                          "Large (10 < X <= 25)", 
                                                          "Medium (5 < X <= 10)", 
                                                          "Small (1 < X <= 5)", 
                                                          "Single (0 < X <= 1)",
                                                          NA))

head(factor(slot(CD8_TCR, "meta.data")$cloneType),20)
saveRDS(CD8_TCR,"cd8_tcr.rds")
unique(CD8_TCR$CTaa)

##delete T cells without TCR data
tcrbarcode=c(combined$HB01_T$barcode,combined$HB03_T$barcode,combined$HB04_T$barcode,combined$HB05_T$barcode,combined$HB06T_T$barcode,
             combined$HB06P_P$barcode,combined$HB07T_T$barcode,combined$HB07P_P$barcode,combined$HB08_T$barcode,combined$HB09T_T$barcode,
             combined$HB09P_P$barcode,combined$HB10T_T$barcode,combined$HB10P_P$barcode,combined$HB11_T$barcode,combined$HB14T_T$barcode,
             combined$HB14P_P$barcode,combined$HB16T_T$barcode,combined$HB16P_P$barcode,combined$HB17T_T$barcode,combined$HB17P_P$barcode)

metadata=CD8_TCR@meta.data[rownames(CD8_TCR@meta.data)%in%tcrbarcode,]

CD8_TCR=CD8_TCR[,rownames(CD8_TCR@meta.data)%in%tcrbarcode]
CD8_TCR@meta.data=metadata

##bar plot  
cell.prop <-as.data.frame(prop.table(table(CD8_TCR$cloneType, CD8_TCR$annotation), margin = 2))
colnames(cell.prop)<-c("celltype","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  geom_bar(stat = "identity",position="fill")+ggtitle("")+theme_classic()+ 
  theme(axis.ticks.length=unit(0.5,'cm'))+guides(fill=guide_legend(title=NULL))+
  RotatedAxis()+scale_fill_d3("category20c")


########BCR
library(patchwork)
add_clonotype2 <- function(tcr_prefix, seurat_obj, type="t"){
  #contig_annotations.csv
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  # Subsets so only the first line of each barcode is kept,as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono)
  #Reorder so barcodes are first column and set them as rownames.
  rownames(tcr) <- tcr[,2]
  tcr[,2] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

#####load data
bcell_bcr<- add_clonotype2("BCR/",bcell_bcr, type="onlyh")

bcell_bcr@meta.data$BCR <- c("No BCR")
bcell_bcr@meta.data$BCR[bcell_bcr@meta.data$b_is_cell==c("TRUE")]<-"Yes"
DimPlot(bcell_bcr, reduction = "umap",group.by = "BCR", 
        repel = TRUE,cols = c("#E6E6E6","#3180D7"))

######
table(bcell_bcr$b_chain)
table(bcell_bcr$b_c_gene)
bcell_bcr@meta.data$BCRtype <- c("nobcr")
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHM")]<-"IgM"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHA1")]<-"IgA1"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHA2")]<-"IgA2"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHD")]<-"IgD"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHE")]<-"IgE"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHG1")]<-"IgG1"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHG2")]<-"IgG2"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHG3")]<-"IgG3"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHG4")]<-"IgG4"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGHM")]<-"IgM"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene==c("IGKC")]<-"kappa"
bcell_bcr@meta.data$BCRtype[bcell_bcr@meta.data$b_c_gene%in%c("IGLC1","IGLC2","IGLC3")]<-"lamda"
DimPlot(bcell_bcr,group.by = "BCRtype")+scale_color_d3("category20c")
table(bcell_bcr$BCRtype)

bcell_bcr@meta.data$igtype <-  bcell_bcr@meta.data$BCRtype
bcell_bcr@meta.data$igtype[bcell_bcr@meta.data$BCRtype%in%c("IgA1","IgA2","IgE",
                                                            "IgG1","IgG2","IgG3","IgG4")]<-"CSIg"
bcell_bcr@meta.data$igtype[bcell_bcr@meta.data$BCRtype%in%c("kappa","lamda","nobcr")]<-"nobcr"
DimPlot(bcell_bcr,group.by = "igtype")
DimPlot(bcell_bcr, reduction = "umap",group.by = "annotation", label = TRUE)+scale_color_d3("category20c")

Idents(bcell_bcr) <- "annotation"
bcell_bcr <- subset(bcell_bcr,idents="pDC",invert=T)
cell.prop <-as.data.frame(prop.table(table(Idents(bcell_bcr), bcell_bcr$igtype), margin = 2))
colnames(cell.prop)<-c("cluster","igtype","proportion")
ggplot(cell.prop,aes(x=factor(cluster,levels = c("Pre-pro_B","Pro_B","Pre_B","Immature_B","Naive_B",
                                                 "Activated_B","Plasma_cell","Plasmablast")),proportion,fill=igtype))+
  geom_bar(stat = "identity",position="fill")+ggtitle("")+theme_bw()+ 
  theme(axis.ticks.length=unit(0.5,'cm'))+guides(fill=guide_legend(title=NULL))+
  RotatedAxis()+scale_fill_d3("category20c")+xlab("")


