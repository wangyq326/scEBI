library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(pheatmap)
#
HB=readRDS("HB.rds")

Idents(HB)="annotation"
marker_gene=FindAllMarkers(HB,only.pos = T)
write.csv(marker_gene,"total_deg.csv")
top10 <- marker_gene %>% group_by(cluster) %>% top_n(10, avg_log2FC)
marker_gene=top10[!duplicated(top10$gene),]
my_order=unique(HB$annotation)

HB2=subset(HB,features=marker_gene$gene)
hb_merge2 <- subset(HB2,
                    subset = annotation %in% my_order)
#
my_order_bk=my_order
#
hb_merge_mat=hb_merge2@assays$RNA@data
hb_merge_mat=as.matrix(hb_merge_mat)
hb_merge_mat[1:4,1:5]
hb_merge_mat=t(hb_merge_mat)
hb_merge_mat=as.data.frame(hb_merge_mat)
hb_merge_mat$annotation=hb_merge2$annotation
#
tmp=aggregate(hb_merge_mat[,-ncol(hb_merge_mat)],by=list(annotation=hb_merge_mat$annotation),FUN=mean)
row.names(tmp)=tmp$annotation
tmp=tmp[my_order,]
tmp=tmp[,-1]
range(tmp)

#order the gene
marker_gene$cluster=factor(marker_gene$cluster,levels = my_order)
marker_gene$cluster
unique(marker_gene$cluster)
marker_gene=marker_gene[marker_gene$cluster%in%my_order,]
marker_gene=marker_gene[order(marker_gene$cluster),]
gene=marker_gene$gene
gene=unique(marker_gene)
#group gene
marker_gene_bk=marker_gene[!duplicated(marker_gene$gene),]
write.csv(marker_gene_bk,'marker_gene_bk.csv',quote = F,row.names = F)
marker_gene_bk_2=read.csv('marker_gene_bk.csv')
#
tmp[1:4,1:5]
tmp_bk=tmp
tmp=tmp[,gene$gene]

#plot
tmp=scale(tmp,center = T,scale = T)
range(tmp)
tmp2=log(tmp+2)
range(tmp2)
tmp2[1:4,1:4]
tmp2=scale(tmp2,center = T,scale = T)
range(tmp2)
#color
#breaks
bk <- c(seq(-3,0,by=0.01),seq(0,8,by=0.01))

tmp2=as.data.frame(t(tmp2))
pheatmap(tmp2,cluster_rows = F,
         color = c(colorRampPalette(colors = c("#32729A","#80B5D6","white"))(length(bk)/4.35),colorRampPalette(colors = c("white","#C00000"))(length(bk)/1.3)),
         cluster_cols = F,cellwidth = 15,cellheight =2,
         show_rownames = F,angle_col = 90,filename = 'my_heat.pdf',height = 40,width = 20
)
