untar("./SNP-alignment-largeDataset+Spar.tgz")
chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI")
load_chromosomes = data.frame() # empty dataframe

for (i in 1:length(chromosomes)) {  #for each chromosome...
  bob = read.csv(paste(
    "./SNP-alignment-largeDataset+Spar/chr" 
    , chromosomes[i],
    "_SNP_scere_largeDataset+Spar.fasta",
    sep = ""), header = FALSE)  #import the file as a data frame,
  StrainID <- grep("^>", bob$V1)  #make a vector of the strain names
  for (j in 1:length(StrainID)) {
    load_chromosomes[j, i] = bob[j*2, 1]  #each strain is a row, each chromosome sequence is a column
  }
}
rownames(load_chromosomes) = bob[StrainID, ]  #name rows with strain names

Csequences <- data.frame(Sequence = rep("", nrow(load_chromosomes)), 
                         row.names = bob[StrainID,])  #new data frame for concatenated sequences
for (k in 1:ncol(load_chromosomes)) {
  Csequences$Sequence <- c(paste(Csequences$Sequence, load_chromosomes[, k]))
}  #concatenate all sequences into a single string

write.csv(Csequences, "./ConcatenatedSequences.csv")  #output