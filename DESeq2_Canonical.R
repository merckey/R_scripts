################################
# 2014-8-13
# Comparrison of EMT-canonical samples : PD2204, PD2412, PD354
# will do calculations for batch control

library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_Canonical_vs_Non-Canonical/Canonical")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_Canonical_vs_Non-Canonical/Canonical/'

sampleFiles <- grep('PD',list.files(directory),value=T)

# view sampleFiles
sampleFiles

# set sampleConditions for Nicole's RNAseq ecad neg/pos experiment
sampleCondition <- c('E_minus','E_plus')
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory = directory, design=~condition)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, levels=c('E_plus','E_minus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)

# since only comparing 1 sample for each condition - multiple hypothesis testing will not work
# instead - I am sorting hits based on lowest p-value
res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# pvalue should be ranked lowest pvalue to high (most sig to least sig)

# save data 'res' to csv!
write.csv(as.data.frame(res),file='2014-8-9-DESeq2-noncanicalwithPD9210.csv')

mcols(res, use.names=T)
write.csv(as.data.frame(mcols(res,use.name=T)),file='2014-8-9-DESeq2-test-conditions-nonCanonical.csv')

ddsMF <- dds
# design(ddsMF) <- formula(~ Batch + condition) # don't need to do this since all batch 1 
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
resMF <- resMF[order(resMF$padj),]

head(resMF)

#order by downregulated genes
resMF.log2 <- resMF[order(resMF$log2FoldChange),]
resMF.down <- resMF.log2
head(resMF.down)
#order by upregulated genes
resMF.up <- resMF[order(-resMF$log2FoldChange),]
head(resMF.up)

resMF.up[1:800,]
write.csv(as.data.frame(resMF.up.214 ),file="2014-9-4-DESeq2_lfc_2fold_up_Canonical.csv")

resMF.down.116 <- resMF.down[1:116,]
write.csv(as.data.frame(resMF.down.411),file="2014-9-4-DESeq2_lfc_2fold_down_Canonical.csv")

write.csv(as.data.frame(resFC),file="2014-9-14-DESeq2_lfc_0.5_Canonical.csv")

resUP <- results(ddsMF, lfcThreshold=0.5, altHypothesis='greater')
resUP <- resUP[order(resUP$log2FoldChange),]
write.csv(as.data.frame(resUP),file="2014-9-14-DESeq2_lfc_up_Canonical.csv")

vsd <- varianceStabilizingTransformation(dds, blind=T)
#heatmap of data
library("RColorBrewer")
library("gplots")
# 200 top fold change genes
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:214]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
par(cex.main=1)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=TRUE, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F,
          main="Canonical EMT transcriptome")
dev.copy(png, "2014-9-5-DESeq2_heatmap_canonical_top221_2.png")
dev.off()

#order by padj
select.padj = order(res$padj,decreasing=F)[1:100]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
par(cex.main=1)
heatmap.2(assay(vsd)[select.padj,], col=my_palette,
          Rowv=F, scale="row", key=TRUE, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F, margins=c(6,6),
          main="Canonical EMT transcriptome")
dev.copy(png, "2014-9-8-DESeq2_heatmap_canonical_padjtop100.png")
dev.off()

# transform raw distrbuted counts for clustering analysis
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# views clustering on individual datasets in unbiased way
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"2014-7-23-DESeq2_variance_stabilizing.png")
# axis is square root of variance over the mean for all samples
dev.off()

# clustering analysis
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-7-23-deseq2_heatmaps.png")
dev.off()


#heatmap of data
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=F)),decreasing=T)[1:3000]
se#hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(assay(rld)[select, ], col=redblue(16),
          Rowv = F, Colv = F, scale= 'non',
          dendrogram = 'none', trace = 'none', margin = c(6,6))
dev.copy(png, '2014-7-23-DESeq2_heatmap10.png')
dev.off()


# PCA analysis
# run principal component analysis on data
# good for visualizing effect of experimental covariats and batch effect
# ideal for examining primary and matching mets
# using plotPCAWithNames.R script - can change number of ntop

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=500)  # can change ntop from 200-500
{
  library("genefilter")
  library("lattice")
  
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))
  
  # extract sample names
  names = colnames(x)
  
  fac = factor(apply( as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  if( nlevels(fac) >= 3 )
    colours = brewer.pal(nlevels(fac), "Dark2")
  else  
    colours = c( "red", "blue" )
  
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    panel=function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
    },
    aspect = "iso", col=colours,
    main = draw.key(key = list(
      rect = list(col = colours),
      text = list(levels(fac)),
      rep = FALSE)))
}

print(plotPCAWithSampleNames(vsd, intgroup=c('condition')))
dev.copy(png, "2014-8-9-DESeq2-PCA-noncanonical.png")
dev.off()

