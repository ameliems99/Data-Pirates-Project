
#Distance Matrix

library(ggplot2)

DataP<-DataDat
DataP$value[DataP$value>0.1]<-NA
ggplot(data = DataP, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+scale_fill_gradientn(colours=c("white","blue","green","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Tree building
DM2 <- DataDM
DM2[DM2>0.2] <- NA
Tree <- njs(DM2)
library(ggtree)
ggtree(Tree, layout = "circular")+geom_tiplab()

#NMDS 

library(vegan)
Dist<- vegdist(DataDM,method="bray",binary=F)
set.seed(20)
NMDSdat<-metaMDS(Dist,k=2,trymax = 100)
DataNM<-data.frame(NMDS1=NMDSdat$points[,1],
                 NMDS2=NMDSdat$points[,2],
                 SampleID=row.names(DataDMmat))
DataNM<-merge(DataNM,Info,by.x="SampleID",by.y="Strain",all.x=T,all.y=F)

library(countrycode)
DataNM$Continent <- countrycode(sourcevar = DataNM[, "Geographic_location.y"],
                            origin = "country.name",
                            destination = "region",
                            nomatch = NULL,
                            custom_match = c('Asia' = 'East Asia & Pacific','Hawaii'='North America','West Africa'='Sub-Saharan Africa'))
qplot(x=NMDS1,y=NMDS2,colour=Continent,data=DataNM,alpha=I(0.4))+theme_bw()


