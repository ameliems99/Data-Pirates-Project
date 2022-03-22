Listed <- list.files("./Scer_SNP_alignment")  #lists all the files in this folder
Keep <- grep("^\\w+_.*.fasta$", Listed)  #vector of the files we want to concatenate
KeepList <- Listed[Keep]  #subset the list of files we want



#load all of the necessary files
FastaList <- list()
library(seqinr)
setwd("./Scer_SNP_alignment")

for (i in 1:length(KeepList)) {
  FastaList[[i]] <- read.fasta(KeepList[i], as.string = TRUE, forceDNAtolower = FALSE)
}

#this loads all of the files into a List of lists (each list item is a list in itself of the sequences for each strain)



setwd("../") #reset working dir to Data-Pirates-Project


