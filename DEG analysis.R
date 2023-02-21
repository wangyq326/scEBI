#####volcano plot
library(ggpubr)
library(ggplot2)
library(ggthemes)

Idents(epc_pure)="sample_type"
res<-FindMarkers(epc_pure)

res$logP<--log10(res$p_val_adj)
head(res)
res$Group="not-significant"
res$Group[which((res$p_val_adj<0.05)&(res$avg_log2FC > 0.8))]="up-regulated"
res$Group[which((res$p_val_adj<0.05)&(res$avg_log2FC < -0.8))]="down-regulated"
table(res$Group)
res$Label<-""
res$Symbol<-rownames(res)
up.genes<-head(res$Symbol[which(res$Group=="up-regulated")],table(res$Group)[[2]])
down.genes<-head(res$Symbol[which(res$Group=="down-regulated")],table(res$Group)[[1]])
deg.top.genes<-c(up.genes,down.genes)
res$Label[match(deg.top.genes,res$Symbol)]<-deg.top.genes

g=ggscatter(res,x="avg_log2FC",y="logP",
            color = "Group",
            title = "EPC:tumor vs peri",
            palette = c("#2f5688","#BBBBBB","#CC0000"),
            size = 2,
            label = res$Label,
            repel = T,
            xlab = "Fold change, log2(gene expression)",
            ylab = "-log10(adj.P value)")+theme_base()+
  geom_vline(xintercept = c(-0.8,0.8),linetype = "dashed")
g

###clusterprofile
library(DOSE)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)

###
data_upregulated<-rownames(res)[res$avg_log2FC > 0]
data_downregulated<-rownames(res)[res$avg_log2FC < 0]

##
data_upregulated<-as.character(data_upregulated)
data_downregulated<-as.character(data_downregulated)

####
test_up<- bitr(data_upregulated,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
head(test_up,2)

test_down<- bitr(data_downregulated,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
head(test_down,2)

##enrich GO 
#ALL
ego_ALL<-enrichGO(gene=test_up$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE) 

barplot(ego_ALL,showCategory = 20,title="EPC HB enriched")

ego_ALL<-enrichGO(gene=test_down$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE) 

barplot(ego_ALL,showCategory = 20,title="EPC liver enriched")

###scMetabolism

library(scMetabolism)
library(ggplot2)
library(rsvd)

countexp.Seurat<-sc.metabolism.Seurat(obj = epc_pure, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score

input.pathway<-(c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", 
                  "Citrate cycle (TCA cycle)",'Pyruvate metabolism',
                  'Pentose phosphate pathway','Fatty acid biosynthesis',
                  'Fatty acid degradation',
                  'Arginine biosynthesis',
                  'Arginine and proline metabolism',
                  "Alanine, aspartate and glutamate metabolism" ,                           
                  "Glycine, serine and threonine metabolism" ,                              
                  "Cysteine and methionine metabolism"))


DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "annotation", norm = "y")
