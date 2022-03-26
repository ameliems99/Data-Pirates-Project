# This code will perform a principal component analysis from the csv file containing the concatenated sequences



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

# In the SNP dataframe, we replace nucleotides indentical to ConsSeq by 0 while the polymorphism are replaced by 1
for (i in 1:147) {
  SNP[(SNP[,i] == ConsSeq),i] = 1
  SNP[(SNP[,i] != 1),i] = 10
}


PCA = princomp(SNP[2:10,2:10])
# meh it's not working

library(ggplot2)


qplot(PCA$scores[,1], PCA$scores[,2])
