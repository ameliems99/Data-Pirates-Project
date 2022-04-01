# 01/04/2022


# This code will perform a principal component analysis from the csv file 
# containing the concatenated sequences

# PCA is a methode widely used in biology: PCA is traditionally used with
# numeric data (eg. gene expression profiles,...) it is used when many 
# parameters where measured for a set of individuals and that it is not possible
# to directly plot the measurements on a 2D graph. 
#
#  
# However it cannot be directly used on DNA sequences. To do so, I chose to 
# convert the DNA sequence in a binary format: Nucleotides will be replaced by a
# 0 if they correspond to the most common nucleotide at this position, and by 1 
# if they differ from it (if they are SNPs) 
# 
# A PCA will be performed on the binary sequences 


#
#-------------------------------------------------------------------------------
library(ggplot2)# import ggplot  because we will need to make graphs

#-------------------------------------------------------------------------------
#The first step is to load the dataset and to split sequences:


#loading of the DataSet, and we keep the length of the concatenated sequences in
# WordLen
DataSeq = read.csv("./ConcatenatedSequences.csv", header = T)
WordLen = nchar(DataSeq[1,2]) # the length of the sequences
# So far the data in DataSeq look like this:


#     1             2
# name_Seq_1 "ACGTACGT"
# name_seq_2 "ACGTACGG"
# name_seq_3 "ACGTACTT"



# creation of a Data frame called SNP, in which each colum contains a sequence 
# split in single nucleotides

SNP = data.frame()

for (i in 1:dim(DataSeq)[1]) {
  SNP[1:WordLen,i] = unlist(strsplit(DataSeq[i,2], ""))
  
}

# now the Data in SNP look like this:

#   1  2   3   4   5   6   7   8
# "A" "C" "G" "T" "A" "C" "G" "T"
# "A" "C" "G" "T" "A" "C" "G" "G"
# "A" "C" "G" "T" "A" "C" "T" "T"

#-------------------------------------------------------------------------------
# Now we want to create a vector which contain the most common nucleotide of 
# each position


# for a given vector x, MostFreq return the most common value of the vector
MostFreq = function(x){
  return(
    names(sort(table(as.character(x)),decreasing=TRUE)[1])
  )
}

# with the MostFreq function we create MFSeq (Sequence w/ most frequent 
# nucleotides): a vector containing the most common nucleotide of each position. 

MFSeq = apply(SNP,1, MostFreq )

# With out example, MFSeq would look like this:

# "A", "C", "G", "T", "A", "C", "G", "T"


#-------------------------------------------------------------------------------
# The goal is now to make a binary version of the SNP dataframe: for each 
#sequence we replace the nucleotide by 0 if it is identical to the most common 
#nucleotide, or by 1 if it is different (whatever the mutation)


# In the BinSNP dataframe, we replace nucleotides indentical to MFSeq 
# by 0 while the polymorphisms are replaced by 1

BinSNP = data.frame()

for (i in 1:147) {
  BinSNP[as.vector(which(SNP[,i] == MFSeq)),i] = 0
  BinSNP[as.vector(which(SNP[,i] != MFSeq)),i] = 1
}

# Now We have a dataframe that look like this 
#
# 0 0 0
# 0 0 0
# 0 0 0
# 0 0 0
# 0 0 0
# 0 0 0
# 0 0 1
# 0 1 0
# 
# because R doesn't like when column are missing, so I inverted the rows and
# columns.
# However, to Perform a PCA we need the the strains as rows and the nucleotides 
# as column so we transpose the data frame. 

PCA_table = t(BinSNP)
#we eventually have a matrix that look like this:

# 0 0 0 0 0 0 0 0 
# 0 0 0 0 0 0 0 1 
# 0 0 0 0 0 0 1 0 


#-------------------------------------------------------------------------------
#We import strain info from the ordered file: the order of the strain in that  
# file is the same as the in sequence file

strain_info = read.csv("ordered_strain_info.csv")

#-------------------------------------------------------------------------------

# we can remove the outgroups and the sequences for which we don't have any info
which(strain_info$Phylogeny == "")
PCA_table = PCA_table[c(1:53,55:146),]
strain_info = strain_info[c(1:53,55:146),]

#-------------------------------------------------------------------------------
# Because we have more nucleotides than samples, it is not possible to use 
# princomp We use prcomp, which is essentially the same thing: the calculation 
# methode is not the same but the result are often very similar

# we calculate the principal components
PCA = prcomp(PCA_table)




# and we plot nice graphs with an appropriate colour code!


colour_code = strain_info$Phylogeny

qplot(PCA$x[,1], PCA$x[,2], colour = colour_code, 
  main = "Genetic distance of yeast strains based on SNPs"
  ) +
  xlab("PC 1") +
  ylab("PC 2") +
  scale_color_discrete("phylogeny") +
  theme_bw()


# Though the data is nicely spread along the PC1, most of the data is on the 
# axis for the PC2, except for West african strains. The way the PC2 is spread 
# is not optimal to read the graph.
# Here I did again the exact same graph but with PC1 and PC3. 


qplot(PCA$x[,1], 
      PCA$x[,3], 
      colour = colour_code, 
      main = "Genetic distance of yeast strains based on SNPs",
      shape =(PCA$x[,2]<=(-20))) +
  
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("phylogeny") +
  scale_shape_discrete("very negative PC2") +
  theme_bw()
#-------------------------------------------------------------------------------

# we can also change the colour code according to substrate of isolation:
colour_code = strain_info$Substrate_of_isolation

qplot(PCA$x[,1],
      PCA$x[,3], 
      colour = colour_code, 
      main = "Genetic distance of yeast strains isolated from different substrate",
      shape =(PCA$x[,2]<=(-20))) +
  
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("substrate of isolation") +
  scale_shape_discrete("very low PC2") +
  theme_bw()

# Or according to the geographic location

colour_code = strain_info$Geographic_location

qplot(PCA$x[,1], 
      PCA$x[,3], 
      colour = colour_code, 
      main = "Genetic distance of yeast strains isolated 
      from different geographic location") + 
  
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("Location") +
  theme_bw()

#-------------------------------------------------------------------------------
# with the squared standard deviations we make a Scree Plot, only with the first 
# 25 bcs otherwise there are too many values and it's unreadable
qplot(x=c(1:25),y=PCA$sdev[1:25]^2, main = "Scree Plot", ylim = c(0, 200)) +
  geom_line() +
  xlab("Component") +
  ylab("Eigenvalue")
