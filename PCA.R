# This code will perform a principal component analysis from the csv file containing the concatenated sequences


# I AM NOT SURE WHAT IS IN THIS CODE IS ACTUALLY VALID, I ALREADY EMAILED MARIA TO ASK HER WHAT SHE THOUGHT ABOUT IT.

#---------------------------------------------------------------------------------------------------------------


Data = read.csv("./ConcatenatedSequences.csv", header = T)

WordLen = nchar(Data[1,2]) # the length of the sequences




SNP = data.frame()

# creation of a Data frame, in which each colum contains a sequence split in single nucleotides
for (i in 1:dim(Data)[1]) {
  SNP[1:WordLen,i] = unlist(strsplit(Data[i,2], ""))
  
}


# for a given vector x, MostFreq return the most common value of the vector
MostFreq = function(x){
  return(
    names(sort(table(as.character(x)),decreasing=TRUE)[1])
  )
}

# with the MostFreq function we create ConsSeq: a vector containing the most common nucleotide of each position. 
ConsSeq = apply(SNP,1, MostFreq )


# In the Binary SNP dataframe, we replace nucleotides indentical to ConsSeq by 0 while the polymorphism are replaced by 1
BinSNP = data.frame()

for (i in 1:147) {
  BinSNP[as.vector(which(SNP[,i] == ConsSeq)),i] = 0
  BinSNP[as.vector(which(SNP[,i] != ConsSeq)),i] = 1
}

# we make the PCA
PCA = princomp(BinSNP, cor = T)

# we pplot the first 2 components
library(ggplot2)
qplot(PCA$scores[1,], PCA$scores[2,])

# TA-DA!!
