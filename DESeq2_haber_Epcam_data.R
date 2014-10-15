# R script for analysis of raw counts from Ting et al 2014, "Singe-Cell RNA sequencing..."
# raw counts downloaded from NCBI gene expression omnibus GSE51372


library("DESeq2")
setwd("~/Desktop/Haber_CTC_data/")
directory <- '/Users/dballi/Desktop/Haber_CTC_data/'
load("2014-10-8-Haber-CTC-WBC-counts.RData")

EpcamCounts <- read.delim("GSE51372_readCounts.Epcam.binned.txt")
EpcamCounts

# set row names with make.names function
row.names(EpcamCounts) <- make.names(EpcamCounts$gene, unique=T)
# remove unnecessary habercounts$symbol coumn
EpcamCounts$gene <- NULL
head(EpcamCounts)

EpcamCounts <- as.matrix(EpcamCounts)
ncol(EpcamCounts)
conds <- factor(c(rep("Epcam+",11), rep("Epcam-",75)))

CountsTable <- data.frame(row.names=colnames(EpcamCounts),conds)

head(CountsTable)

dds <- DESeqDataSetFromMatrix(countData=EpcamCounts, 
                              colData=CountsTable, 
                              design=~conds)

colData(dds)$conds <- factor(colData(dds)$conds, 
                                      levels=c('Epcam+','Epcam-'))


dds <- DESeq(dds)

res <- results(dds)
table(res$padj < 0.5)
res <- res[order(res$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(assay(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="2014-10-8-results-Ha-DESeq2.csv")
reslfc <- subset(resdata, resdata$padj < 0.05)
reslfc <- reslfc[order(-reslfc$log2FoldChange),]
write.csv(reslfc, file="2014-10-10-lfc-epcam-DESeq2.csv")

resabslfc <- res[order(-res$log2FoldChange),]
write.csv(resabslfc, file="2014-10-10-absolute-lfc-epcam-DESeq2.csv")

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
text3d(pc$x, text=colnames(EpcamCounts), adj=-0.2, font=1, cex=0.8, col=cond_color)
# save as postscript pdf file then edit in gimp
rgl.snapshot("2014-10-8-Haber-Epcam+_Epcam-_3Dpca.png",'png')


save.image("~/Desktop/Haber_CTC_data/2014-10-10-Epcam_comparisions.RData")


