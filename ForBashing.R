untar("./Scer_SNP_alignment.tgz")

ChrI <- read.csv("./Scer_SNP_alignment/ChrI_Scer_SNP_alignment.fasta", header = FALSE) #read as data frame

StrainID <- grepl("^>", ChrI$V1) #make a vector of the strain names
DF<- data.frame(StrainID = ChrI[StrainID == TRUE,],  #make a separate column for strain names
                Sequence = ChrI[StrainID == FALSE,]) #column for just the sequences 

View(DF)


##########

chromosomes = c("I", "II")

all_chromosomes_set = data.frame() # empty dataframe
library(seqinr)
for (i in 1:length(chromosomes)) { # for each chromosome (using length is more "reproducible")
  #we read the fasta file 
  a_chromosome = read.fasta(paste(
    "./SNP-alignment-largeDataset+Spar/chr"
    , chromosomes[i],
    "_SNP_scere_largeDataset+Spar.fasta",
    sep = ""), as.string = TRUE, forceDNAtolower = F)
  for (j in 1:10){ # and for the 10 first sequence we add the concatenated sequence as a string in the dataset.
    all_chromosomes_set[1, j] = paste(all_chromosomes_set[1, j], unlist(a_chromosome[j]), sep = "", collapse = "")  
    #when set to [1, j] instead of [i, j], the whole sequence gets pasted inthe same cell (instead of different rows for each chromosome)
  }
}
colnames(all_chromosomes_set) = attr(a_chromosome,"name")[1:10]
rownames(all_chromosomes_set) = "Sequences"

View(all_chromosomes_set)
