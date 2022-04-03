# The purpose of this code is to output a nice csv file with 
#the name of strain, their position in the concatenated sequence.csv

library(seqinr) # allows to access read.fasta // you probably have to install it first

# we read a short csv file and we get the names of the strains for which we have SNPs 
untar("./SNP-alignment-largeDataset+Spar.tgz")
Fastfile = read.fasta("./SNP-alignment-largeDataset+Spar/chrI_SNP_scere_largeDataset+Spar.fasta")
all_strains = names(Fastfile)

# We get load the file containing the info about SNPs
Data = read.csv("strain_info.csv")

# we create a new data frame
Info = data.frame()
for ( i in 1:length(all_strains)) { # for every strain in the SNPs files...
  Info[i, 1:8] = Data[(Data$Strain == all_strains[i]), ] #we add a row t Data2 containing the corresponding info for the strain
}

#Regroup substrate of isolation into categories
Info$Substrate_of_isolation <- gsub(".*(Wine|wine).*", "Wine", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub("(Quercus|Quecus|Oak).*", "Oak", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(Fermentation|fermentation|Fermented|Brewery|brewing|
                                    brewerie|beer|sake|Bioethanol|bioethanol|Champagne|ethanol).*", 
                                    "Other fermentation", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(Fruit|fruit|Grape|fig|Apple|apple|Vineyard).*", "Fruit", 
                                    Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(palm|Ficus|Opuntia|Fraxinus|Fagus|Castanea).*", "Other tree",
                                    Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(Sugar|Coconut).*", "Other plant material", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*[^Wine|Oak|Other fermentation|Fruit|Other tree|
                                    Other plant material].*", "Other", Info$Substrate_of_isolation)

#Regroup geographic location by country 
Info$Geographic_location <- gsub(".*,\\s(\\w+$)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(Greece)", "\\1", Info$Geographic_location)
Info$Geographic_location <- gsub(".*(Portugal).*", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(Greece)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(South Africa)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(UK)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(France)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub("(Asia).*", "\\1", Info$Geographic_location)


# We output the data in the following file:
write.csv(Info, "ordered_strain_info.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------------

#This script reformats data from the chromosome I fasta file into .csv format for analysis

load_chromosome = data.frame() # empty dataframe

bob1 <- read.csv("./SNP-alignment-largeDataset+Spar/chrI_SNP_scere_largeDataset+Spar.fasta", header = FALSE)  #import the file as a data frame,
StrainID <- grep("^>", bob1$V1)  #make a vector of the strain names
for (j in 1:length(StrainID)) {
   load_chromosome[j, 1] = bob1[j*2, 1]  #each strain is a row, each chromosome sequence is a column
}
rownames(load_chromosome) = bob1[StrainID, ]  #name rows with strain names
colnames(load_chromosome) <- "Sequence"

write.csv(load_chromosome, "./ChrISequences.csv")  #output

#----------------------------------------------------------------------------------------------

#This part of the script generates a distance matrix, phylogeny, and NMDS


#Load requires libraries
library(BiocManager)
library(Biostrings)
library(muscle)
library(ape)
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggtree)

#load data
MyData<-read.csv("./ChrISequences.csv")

#Clean up Data
MyData$Sequence<-gsub("[^ATCG]", "", MyData$Sequence)

#Removing S. paradoxus species 
MyData2 <- MyData[!grepl("Spar", MyData$X),]

#Aligning sequences with new data not containing the S. paradoxus 
DataDNAstring<-DNAStringSet(MyData2$Sequence)
MyData2$X <- gsub(">", "", MyData2$X)
names(DataDNAstring)<-paste(MyData2$X)
DataAlign<-muscle::muscle(stringset = DataDNAstring, quiet = T)


#Distance Matrix
DataAlignBin<-as.DNAbin(DataAlign)
DataDM<-dist.dna(DataAlignBin, model = "K80")
DataDMmat<-as.matrix(DataDM)
DataDat<-melt(DataDMmat)

#Distance Matrix Plot
png("./images/DistanceMatrix.png")
  ggplot(data = DataDat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()+scale_fill_gradientn(colours=c("white","blue","green","red")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("Yeast Strains") + 
    ylab("Yeast Strains") + 
    labs(title = "Distance Matrix of Yeast Strains")
dev.off()


#Phylogeny plot
Tree <- nj(DataDM)
Info <- read.csv("./ordered_strain_info.csv") #import info data
Info <- Info[!Info$Strain == "Spar",] #remove S. paradoxus

png("./images/Phylogeny.png")
  ggtree(Tree, branch.length = 'none', layout = "circular") %<+% Info + 
    geom_tree(aes(colour = Substrate_of_isolation)) + 
    labs(title = "Phylogeny of Yeast Strains", colour = "Substrate of Isolation") 
dev.off()


#NMDS plot

library(vegan)
library(countrycode)

Dist<- vegdist(DataDM, method="bray", binary=FALSE)
set.seed(20)
NMDSdat<-metaMDS(Dist,k=2,trymax = 100)
DataNM<-data.frame(NMDS1=NMDSdat$points[,1],
                   NMDS2=NMDSdat$points[,2],
                   SampleID=row.names(DataDMmat))
DataNM<-merge(DataNM,Info,by.x="SampleID",by.y="Strain",all.x=T,all.y=F)

DataNM$Continent <- countrycode(sourcevar = DataNM[, "Geographic_location"],
                                origin = "country.name",
                                destination = "region",
                                nomatch = NULL,
                                custom_match = c('Asia' = 'East Asia & Pacific','Hawaii'='North America','West Africa'='Sub-Saharan Africa'))

png("./images/NMDS_Continent.png")
  qplot(x=NMDS1,y=NMDS2,colour=Continent,data=DataNM,alpha=I(0.4))+theme_bw()
dev.off()

png("./images/NMDS_substrate.png")
  qplot(x=NMDS1,y=NMDS2,colour=Substrate_of_isolation,data=DataNM,alpha=I(0.4))+theme_bw()
dev.off()