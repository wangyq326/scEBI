merge_meta = HB@meta.data
cell_type = levels(factor(merge_meta$annotation))
cell_type

a <- merge_meta %>% group_by(annotation,sample_type) %>%summarise(a=n()) 
ac <- merge_meta %>% group_by(sample_type) %>% summarise(ac=n())  
ab <- merge_meta %>% group_by(annotation) %>% summarise(ab=n()) 
abcd = dim(merge_meta)[1]  
dist.data <- merge(a, ac)
dist.data <- merge(dist.data, ab)

dist.data['b'] = dist.data$ab - dist.data$a       
dist.data['c'] = dist.data$ac - dist.data$a
dist.data['d'] = abcd - dist.data$ab - dist.data$c

head(dist.data)

dist.data[,c('p.value','Ea', 'Ec', 'Eb', 'Ed')] = t(apply(dist.data, 1, function(x){
  x=as.numeric(x[c('a', 'c', 'b', 'd')])
  k.test = chisq.test(matrix(x, nrow=2,ncol=2))
  return(c(k.test$p.value, as.vector(k.test$expected)))
}))
dist.data['Reo'] = dist.data$a / dist.data$Ea

library(ggplot2)
library(reshape2)
library(pheatmap)
dist.heatmap.data <- dist.data[,c("annotation", 'sample_type', 'Reo','p.value')]

dist.heatmap<-data.frame(matrix(data=0, nrow = length(cell_type), ncol=3))
rownames(dist.heatmap) = cell_type
colnames(dist.heatmap) = c('tumor',"non_tumor",'p.value')
for(i in 1:dim(dist.heatmap.data)[1]){
  ct = as.character(dist.heatmap.data[i, 'annotation'])
  ts = dist.heatmap.data[i, 'sample_type']
  reo = dist.heatmap.data[i, 'Reo']
  p=dist.heatmap.data[i, 'p.value']
  dist.heatmap[ct, ts] = reo
  dist.heatmap[ct, 3] = p}


dist.heatmap1<-dist.heatmap[,1:2]

rownames(dist.heatmap1)
annote.heatmap = dist.heatmap
annote.heatmap$annote[annote.heatmap$p.value>=0.01] = ''
annote.heatmap$annote[annote.heatmap$p.value<0.01] = '*'
annote.heatmap$annote[annote.heatmap$p.value<0.001] = '**'
annote.heatmap$annote[annote.heatmap$p.value<0.0001] = '***'
annote.heatmap<-annote.heatmap %>% arrange(tumor) %>% arrange(non_tumor)

annote.heatmap$tumor[annote.heatmap$tumor<1] = ""
annote.heatmap$tumor[annote.heatmap$tumor>1] = annote.heatmap$annote[annote.heatmap$tumor>1]
annote.heatmap$non_tumor[annote.heatmap$non_tumor<1] = ""
annote.heatmap$non_tumor[annote.heatmap$non_tumor>1] = annote.heatmap$annote[annote.heatmap$non_tumor>1]

annote.heatmap<-annote.heatmap[,1:2]

rownames(dist.heatmap1)<-factor(rownames(dist.heatmap1),level=cell_type)

dist.heatmap1<-dist.heatmap1 %>% arrange(tumor) %>% arrange(non_tumor)

cell_cluster_dist=pheatmap(dist.heatmap1,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           treeheight_row = 0,
                           color = colorRampPalette(c("#FFFFFF","#E73f2e"))(50),
                           display_numbers =as.matrix(annote.heatmap),
                           #display_numbers =round(as.matrix(dist.heatmap1), 2),
                           cellwidth = 30,cellheight = 16, fontsize = 10,
                           border_color = '#ffffff',
                           angle_col='45',
                           main = 'cell type',
                           breaks = seq(0,2.5,by=2.5/50)
)