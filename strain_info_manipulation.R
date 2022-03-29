# The purpose of this code is to output a nice csv file with 
#the name of strain, their position in the concatenated sequence.csv

#-------------------------------------------------------------------------------



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
  Info[i,1:8] = Data[(Data$Strain == all_strains[i]),] #we add a row t Data2 containing the corresponding info for the strain
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
Info$Substrate_of_isolation <- gsub(".*[^Wine|Oak|Other fermentation|Fruit|Other tree|Other plant material].*",
                                    "Other", Info$Substrate_of_isolation)



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
write.csv(Info, "ordered_strain_info.csv",row.names = F)
