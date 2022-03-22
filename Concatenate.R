Listed <- list.files("./Scer_SNP_alignment")
Keep <- grep("^\\w+_.*.fasta$", Listed)
KeepList <- Listed[Keep]

FastaList <- list()
library(seqinr)
setwd("./Scer_SNP_alignment")

for (i in 1:length(KeepList)) {
  FastaList[[i]] <- read.fasta(KeepList[i], as.string = TRUE, forceDNAtolower = FALSE)
}



setwd("../")