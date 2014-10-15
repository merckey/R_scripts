# R script for analysis of raw counts from Ting et al 2014, "Singe-Cell RNA sequencing..."
# raw counts downloaded from NCBI gene expression omnibus GSE51372


library("DESeq2")
setwd("~/Desktop/Haber_CTC_data//")
directory <- '/Users/dballi/Desktop/Haber_CTC_data/'
load("2014-10-8-Haber-CTC-WBC-counts.RData")

GSE51372_readCounts.simple <- read.delim("GSE51372_readCounts.simple.txt", row.names='symbol')
HaberCounts <- GSE51372_readCounts.simple

# set row names with make.names function
row.names(HaberCounts) <- make.names(HaberCounts$symbol, unique=T)
# remove unnecessary habercounts$symbol coumn
HaberCounts$symbol <- NULL
head(HaberCounts)

# sub selecting on WBC, and single cell CTCS "MP named"
x <- HaberCounts[,grep("WBC",conames(HaberCounts),value=T)]
y <- HaberCounts[,grep("MP",colnames(HaberCounts),value=T)]
xy <- data.frame(x,y)

head(xy)
# remove primary Tumor samples
xy <-xy[,-grep("Tu",colnames(xy))]

# rename to HaberDF
HaberDF <- xy

write.csv(HaberDF, file="2014-10-8-Haber-CTC-WBC.csv")
save.image("~/Desktop/Haber_CTC_data/2014-10-8-Haber-CTC-WBC-counts.RData")

# so - we have a data frame of each single cell sequenced 
# for WBC and CTC (lineage labeled and non-lineage labelled)
# for each gene in the genome 

head(HaberDF)

countdata <- HaberDF

countdata <- as.matrix(countdata)

conds <- factor(c(rep("WBC",12), rep("CTC",93)))

coldata <- data.frame(row.names=colnames(countdata),conds)

head(coldata)

dds <- DESeqDataSetFromMatrix(countData=countdata, 
                              colData=coldata, 
                              design=~conds)

colData(dds)$conds <- factor(colData(dds)$conds, 
                                      levels=c('WBC','CTC'))


dds <- DESeq(dds)

res <- results(dds)
table(res$padj < 0.5)
res <- res[order(res$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(assay(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="2014-10-8-results-Haber-DESeq2.csv")
reslfc <- subset(resdata, resdata$padj < 0.05)
reslfc <- reslfc[order(-reslfc$log2FoldChange),]
write.csv(reslfc, file="2014-10-8-lfc-Haber-DESeq2.csv")

resabslfc <- res[order(-res$log2FoldChange),]
write.csv(resabslfc, file="2014-10-8-absolute-lfc-Haber-DESeq2.csv")

vsd <- varianceStabilizingTransformation(dds, blind=T)

# 3D PCA plot and 2D PCA
library(grDevices)
library(rgl)
library("genefilter")
library("ggplot2")

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

scores <- data.frame(pc$x, conds, stringsAsFactors=T)
head(scores)
pca <- ggplot(scores, aes(x=PC1, y=PC2,col=(factor(conds))))  # or factor(grouping)
(pca <- pca + geom_point(size = 5))
(pca <- pca + ggtitle("Principal Components"))
(pca <- pca + scale_colour_brewer(name = " ", palette = "Set1"))
(pca <- pca + theme(
  plot.title = element_text(face='bold'),
  legend.justification=c(0,0),
  legend.position=c(0,0), 
  legend.key = element_rect(fill='NA'),
  legend.text = element_text(size=10, face="bold"),
  axis.text.y=element_text(colour="Black"), 
  axis.text.x=element_text(colour="Black"),
  axis.title.x=element_text(face="bold"), 
  axis.title.y=element_text(face='bold'),
  panel.grid.major.x=element_blank(), 
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x=element_blank(), 
  panel.grid.minor.y=element_blank(), 
  panel.background=element_rect(color='black',fill=NA)
))

# 3D pca
# 3D pca with rgl plot3d using same pc 
plot3d(pc$x,size=10,lit=T, box=F, alpha=0.75, axes=T)
cond_color <- as.integer(conds)
text3d(pc$x, text=colnames(countdata), adj=-0.2, font=1, cex=0.8, col=cond_color)
# save as postscript pdf file then edit in gimp
rgl.snapshot("2014-10-8-Haber-CTCvsWBC-3Dpca.png",'png')


save.image("~/Desktop/Haber_CTC_data/2014-10-8-Haber-CTC-WBC-counts.RData")


