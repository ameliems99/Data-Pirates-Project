# The purpose of this code is to output a nice csv file with 
#the name of strain, their position in the concatenated sequence.csv

#-------------------------------------------------------------------------------



library(seqinr) # allows to access read.fasta // you probably have to install it first


# we read a short csv file and we get the names of the strains for which we have SNPs 
Fastfile = read.fasta("doi_10.5061_dryad.hm2jf__v1/SNP-alignment-largeDataset+Spar/SNP-alignment-largeDataset+Spar/chrI_SNP_scere_largeDataset+Spar.fasta")
all_strains = names(Fastfile)


# We get load the file containing the info about SNPs
Data = read.csv("strain_info.csv")

# we create a new data frame
Data2 = data.frame()


for ( i in 1:length(all_strains)) { # for every strain in the SNPs files...
  Data2[i,1:8] = Data[(Data$Strain == all_strains[i]),] #we add a row t Data2 containing the corresponding info for the strain
}


# We output the data in the following file:
write.csv(Data2, "ordered_strain_info.csv",row.names = F)
