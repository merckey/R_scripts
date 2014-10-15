# using EBseq on RSEM counts data 

library(EBSeq)
setwd(dir = "Desktop/RNAseq_Nicole_Ecad/DATA_RSEM/")
directory <- "~/Desktop/RNAseq_Nicole_Ecad//DATA_RSEM"

df <- read.delim(file = "2014-10-14-RSEM_ExpectedCounts_genes.txt")
row.names(df) <- make.names(df$X, unique=T)
df$X <- NULL
head(df)

df <- as.matrix(df)
df <- round(df) # as RSEM calculates isoforms are partial counts 
# need hwole numbers for DEseq
conds <- c("Ecad-","Ecad+","Ecad-","Ecad+",
           "Ecad-","Ecad+","Ecad-","Ecad+",
           "Ecad+","Ecad-","Ecad-","Ecad+",
           "Ecad-","Ecad+","Ecad-","Ecad+")

Sizes=MedianNorm(df)

EBout <- EBTest(Data=df, Conditions = conds, sizeFactors=Sizes,maxround=5)

PP <- GetPPMat(EBout)
str(PP)
head(PP)

DEfound_5 = rownames(PP)[which(PP[,"PPDE"]>=0.95)]
str(DEfound)
head(DEfound)

# calculate fold change
GeneFC <- PostFC(EBout)
str(GeneFC)
PlotPostVsRawFC(EBout,GeneFC)
c