Idents(HB)='annotation'
epc_pure=subset(HB,idents = c("Late_Erythroid","Mid_Erythroid","Early_Erythroid"))

##cell phse
epc_pure <-CellCycleScoring(epc_pure,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
DimPlot(epc_pure,group.by = 'Phase',split.by = "sample_type")
ggplot(epc_pure@meta.data,aes(sample_type,S.Score,fill=sample_type))+
  geom_boxplot()+theme_classic()+stat_compare_means()+NoLegend()


##cytotrace
epc_pure$group=paste0(epc_pure$annotation,epc_pure$sample_type)
phe=as.character(epc_pure$group)
names(phe)=rownames(epc_pure@meta.data)
mat_tumor<- as.data.frame(epc_pure@assays$RNA@counts)
mat_tumor[1:4,1:4]
results <- CytoTRACE(mat = mat_tumor,subsamplesize = 1000,ncores = 10)


###reference: https://github.com/seqyuan/CytoTRACE