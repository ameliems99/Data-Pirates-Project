
#Distance Matrix

library(ggplot2)

PDat2<-DataDat
PDat2$value[PDat2$value>0.1]<-NA
ggplot(data = PDat2, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+scale_fill_gradientn(colours=c("white","blue","green","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Tree building
DM2 <- DataDM
DM2[DM2>0.2] <- NA
Tree <- njs(DM2)
library(ggtree)
ggtree(Tree, layout = "circular")+geom_tiplab()




OTU_dist<- vegdist(OTUs,method="bray",binary=F)
OTU_tree2<-nj(OTU_dist)
ggtree(OTU_tree2,layout="rectangular") %<+% Samples + 
  geom_tiplab(aes(colour=Species)) + 
  theme(legend.position="right")


OTU_bin_dist<-dist(DataAlignBin,method='K80')
library(ape)
library(ggtree)
OTU_tree<-nj(OTU_bin_dist)
ggtree(OTU_tree,layout="rectangular")


library(vegan)
set.seed(13)
NMDSdat<-metaMDS(DataDM,k=2,trymax = 1000)

PDat<-data.frame(NMDS1=NMDSdat$points[,1],
                 NMDS2=NMDSdat$points[,2],
                 SampleID=row.names(DataDM))
PDat<-merge(PDat,Info,by="SampleID",all.x=T,all.y=F)  
qplot(x=NMDS1,y=NMDS2,colour=Substrate_of _isolation,alpha=I(0.6),data=PDat)+theme_bw()
