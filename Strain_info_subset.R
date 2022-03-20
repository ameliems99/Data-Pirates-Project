ChrI <- read.csv("./Scer_SNP_alignment/chrI_Scer_SNP_alignment.fasta", header = FALSE)
StrainID <- grepl("^>", ChrI$V1)
DF<- data.frame(StrainID = ChrI[StrainID == TRUE,], 
                        Sequence = ChrI[StrainID == FALSE,])
DF$StrainID <- gsub(">", "", DF$StrainID)

Strain_info <- read.csv("./Strain_info.csv", header = TRUE)
This_study <- Strain_info[Strain_info$Reference == "This study",]

Vec <- DF$StrainID %in% This_study$Strain
Vec2 <- DF$StrainID %in% This_study$Other_designations

DF_TS <- DF[Vec == TRUE,]
DF_OD <- DF[Vec2 == TRUE,]
DF3 <- rbind(DF_TS, DF_OD)