# heatmap for Gbx2 high group
# 2014-8-22

library("DESeq2")
library("gplots")
load("/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_GBX2_group/.RData")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_GBX2_group/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTseq_DESeq2_analysis_GBX2_group'

# read table of significantly expressed gene raw counts

data <- read.table("sig.counts.txt",header=T,row.names='gene', stringsAsFactors=F)
data2 <- read.table('sig.counts2.txt',header=T,row.names='gene',stringsAsFactors=F)
data <- na.omit(data)
data2<- na.omit(data2)
# remove any rows with 0 raw reads 
data <- data[apply(data!=0,1,any), , drop=F]
data2 <- data2[apply(data!=0,1,any), , drop=F]

dist2 <- dist(data,method='maximum',diag=F,upper=F,p=2)

heatmap(as.matrix(data2), cexCol=0.75,labRow=NA)

heatmap.2(assay(vsd)[select,], col=redgreen(75),  scale="row", key=T, keysize=1,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)
dev.copy(png,"2014-8-22-heatmap2-Gbx2group.png")
dev.off()

data2 <- log2(data2)
boxplot(data,cex.axis=0.8,las=2,main='Original distrubtion of data',ylab="Log2(intensity")

library(preprocessCore)
data2.nm <- normalize.quantiles(as.matrix(data2))
data.nm <- na.omit(data.nm)

rownames(data2.nm) = rownames(data2)
colnames(data2.nm) = colnames(data2)

boxplot(data.nm)

library(limma)
design = cbind(Eplus = c(1,1,1,0,0,0),
              Eminus = c(0,0,0,1,1,1))

fit = lmFit(data, design=design) # not normalized data for now

contrastMatrix = makeContrasts("Eplus","Eminus", levels=design)

fit2 <- contrasts.fit(fit,contrasts=contrastMatrix)
fit2 <- eBayes(fit2)

data3 <- as.matrix(fit) # may need to na.omit
data4 <- as.matrix(fit2) # may need to na.omit
Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
          rep("brown",323))

heatmap.2(assay(vsd)[select,], col=redblue(256), dendrogram="both",
          scale="row", key=T, keysize=0.5, density.info="none",
          trace="none",cexCol=1.2, labRow=NA, RowSideColors=Label,
          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
          lwid=c(1.5,0.2,2.5,2.5))
