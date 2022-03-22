install.packages("seqinr")
library(seqinr)
ChrI <- read.fasta("./Scer_SNP_alignment/chrI_Scer_SNP_Alignment.fasta", as.string = TRUE, forceDNAtolower = FALSE)