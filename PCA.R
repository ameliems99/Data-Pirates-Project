# This code will perform a principal component analysis from the csv file containing the concatenated sequences



library(ggplot2)# import ggplot to make nice graphs

#---------------------------------------------------------------------------------------------------------------
#The first step is to load the dataset and to format it:


#loading of the DataSet, and we keep the length of the concatenated sequences in
# WordLen
Data = read.csv("./ConcatenatedSequences.csv", header = T)
WordLen = nchar(Data[1,2]) # the length of the sequences
# So far the data in Data look like this:


#     1             2
# name_Seq_1 "ACGTACGT"
# name_seq_2 "ACGTACGG"
# name_seq_3 "ACGTACTT"



# creation of a Data frame called SNP, in which each colum contains a sequence split in single nucleotides
SNP = data.frame()

for (i in 1:dim(Data)[1]) {
  SNP[1:WordLen,i] = unlist(strsplit(Data[i,2], ""))
  
}

# now the Data in SNP look like this:

#   1  2   3   4   5   6   7   8
# "A" "C" "G" "T" "A" "C" "G" "T"
# "A" "C" "G" "T" "A" "C" "G" "G"
# "A" "C" "G" "T" "A" "C" "T" "T"

#-------------------------------------------------------------------------------
# Now we want to create a vector which contain the most common nucleotide of each 
#position


# for a given vector x, MostFreq return the most common value of the vector
MostFreq = function(x){
  return(
    names(sort(table(as.character(x)),decreasing=TRUE)[1])
  )
}

# with the MostFreq function we create ConsSeq: a vector containing the most common nucleotide of each position. 
ConsSeq = apply(SNP,1, MostFreq )

# With out example, ConsSeq would look like this:

# "A", "C", "G", "T", "A", "C", "G", "T"


#-----------------------------------------------------------------------------------
# The goal is now to make a binary version of the SNP dataframe: for each sequence we 
# replace the nucleotide by 0 if it is identical to the most common nucleotide, or by 1 
# if it is different (whatever the mutation)


# In the BinSNP dataframe, we replace nucleotides indentical to ConsSeq 
# by 0 while the polymorphisms are replaced by 1

BinSNP = data.frame()

for (i in 1:147) {
  BinSNP[as.vector(which(SNP[,i] == ConsSeq)),i] = 0
  BinSNP[as.vector(which(SNP[,i] != ConsSeq)),i] = 1
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
# To Perform a PCA, we need the the strains as rows and the nucleotides as column
# so we transpose the data frame. 

PCA_table = t(BinSNP)
#we eventually have a matrix that look like this:

# 0 0 0 0 0 0 0 0 
# 0 0 0 0 0 0 0 1 
# 0 0 0 0 0 0 1 0 


#----------------------------------------------------------------------------------
#We import strain info from the ordered file: the order of the strain in that file 
# is the same as the in sequence file

strain_info = read.csv("ordered_strain_info.csv")

#-----------------------------------------------------------------------------------

# we can remove the outgroups and the sequences for which we don't have any info.
which(strain_info$Phylogeny == "")
PCA_table = PCA_table[c(1:53,55:146),]
strain_info = strain_info[c(1:53,55:146),]

#-----------------------------------------------------------------------------------
# Because we have more nucleotides than samples, it is not possible to use princomp
# We use prcomp, which is essentially the same thing: the calculation methode is not the same
# but the result are often very similar

# we calculate the principal components
PCA = prcomp(PCA_table)


# and we plot nice graphs with an appropriate colour code!
colour_code = strain_info$Phylogeny

qplot(PCA$x[,1], PCA$x[,2], colour = colour_code, main = "PCA of the SNP") +
  xlab("PC 1") +
  ylab("PC 2") +
  scale_color_discrete("phylogeny") +
  theme_bw()


# Though the data is nicely spread along the PC1, most of the data is on the axis for the 
# PC2, except for West african strains. The way the PC2 is spread is not optimal to read the
# graph.
# Here I did again the exact same graph but with PC1 and PC3. 


qplot(PCA$x[,1], PCA$x[,3], colour = colour_code, main = "PCA of the SNP") +
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("phylogeny") +
  theme_bw()
#-------------------------------------------------------------------------------

# we can also change the colour code according to substrate of isolation:
colour_code = strain_info$Substrate_of_isolation

qplot(PCA$x[,1], PCA$x[,3], colour = colour_code, main = "PCA of the SNP") +
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("phylogeny") +
  theme_bw()

# Or according to the geographic location

colour_code = strain_info$Geographic_location

qplot(PCA$x[,1], PCA$x[,3], colour = colour_code, main = "PCA of the SNP") +
  xlab("PC 1") +
  ylab("PC 3") +
  scale_color_discrete("phylogeny") +
  theme_bw()

#-------------------------------------------------------------------------------
# with the squared standard deviations we make a Scree Plot, only with the first 10 bcs otherwise there are too
#many values and it's unreadable
qplot(x=c(1:10),y=PCA$sdev[1:10]^2) +
  geom_line() +
  xlab("Component") +
  ylab("Eigenvalue")
