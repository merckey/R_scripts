### bioinformatic analysis of RNA-seq data using DESeq2
### recipe is from Dave Wheeler's blog at Massey University
### http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###
### analysis per DATE on htseq count -intersections-nonempty
### controlling for batch
### removing PD8640 as that sample is really weird 
### also removing PD1849
### excluding PD3044 as it does not look like good sort

### RNAseq batch 1, 2 and 3 for EMT RNAsequencing experiment 
### including new batch on 2015-2-11

library("DESeq2")
library("genefilter")
library("reshape2")
library("ggplot2")
library('heatmap3')
library("RColorBrewer")

setwd("~/Documents/RNAseq_Nicole_Ecad/DESeq2_analysis_batch3/")
directory <- '/Users/dballi/Documents/RNAseq_Nicole_Ecad/Counts/'
save.image("2015-2-11-Batch3-DESeq2.RData")
load("2015-2-11-Batch3-DESeq2.RData")
# can merge individual sample files
sampleFiles <- grep('PD',list.files(directory),value=T)

# view sampleFiles
sampleFiles

sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch1","Batch1","Batch3","Batch3","Batch3","Batch3","Batch3","Batch3",
                 "Batch1","Batch1","Batch1", "Batch1",'Batch2','Batch2')


# set sampleConditions for Nicole's RNAseq ecad neg/pos  comparisions
cdh1_status <- c('E_minus','E_plus')


# set LRT conditions for three way comparision
# L1 - Canonical - Eminus
# L2 - Noncanonical - Eminus
# L3 - Eplus

#conds <- c("L1","L3",
 #          "L1","L3",
  #         "L2",'L3',
   #        "L2",'L3',
  #         "L1","L3",
  #         "L2",'L3',
  #         "L1",'L3',
  #         "L2",'L3',
  #         "L2",'L3',
  #         "LX",'LX')

# also ID sample conditions 
sampleCons <- c("Canonical", 'Canonical', 'Non-Canonical', 'Non-Canonical',
                'Non-Canonical', 'Non-Canonical', 'Canonical', 'Canonical',
                'Non-Canonical', 'Non-Canonical', 'Non-Canonical', 'Non-Canonical',
                'Non-Canonical', 'Non-Canonical', 'Non-Canonical', 'Non-Canonical',
                'Canonical', 'Canonical', 'Non-Canonical', 'Non-Canonical', 
                'Non-Canonical', 'Non-Canonical')



# set sample table with each samples + condition + cdh1_status
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          samplebatch = sampleBatch,
                         # condition = conds,
                          samplecons = sampleCons,
                          cdh1_status = cdh1_status)

# view sampleTable
sampleTable 

################################################################################
#
# 
# Perform DE analysis on E vs M for full group - not controling for batch
#
#
################################################################################

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                       directory = directory, 
                                       design= ~cdh1_status) 

colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$cdh1_status, 
                                      levels=c('E_plus','E_minus'))

ddsEvsM <- DESeq(ddsHTseq)

# controlling for batch effect
design(ddsEvsM) <- formula(~ samplebatch + condition)
ddsEvsM <- DESeq(ddsEvsM)
resEvsM <- results(ddsEvsM)
resEvsM <- resEvsM[order(resEvsM$padj),]
head(resEvsM)
table(resEvsM$padj < 0.1)
resdata <- merge(as.data.frame(resEvsM), as.data.frame(counts(ddsEvsM,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'Gene'
head(resdata)
write.csv(resdata, file="2015-2-11-DESeq2-batch3-batchcontrol-results-with-normalied.csv")

# variance stabilization normalization
vsd <- varianceStabilizingTransformation(ddsEvsM,blind=T)
write.csv(assay(vsd), file = "2015-2-11-DESeq2-batch3-variancestabilized.csv")

################################################################################
#
# 
# EDA of E vs M group
#
#
################################################################################

# maplot 

dds <- ddsEvsM
plotMA(dds, ylim=c(-10,10),main = "Ecad+/- RNAseq batch affect controlled")
dev.copy(pdf, "2015-2-11-DESeq2_MAplot_batch3.pdf")
dev.off()


# Unsupervised heirarchical clustering with complete linkage and 1 minus the Pearson's correlation coefficient 
# complete linkage 
hclust_comp <- function(x,method='complete')
  hclust(x,method=method)

# Pearson's correlation coefficient for distance metric
dist_comp <- function(x)
  as.dist(1-cor(t(x)))

# HEATMAP
# add more for new samples
colsidecolors = c("#D95F02","#D95F02",
                  "#1B9E77","#1B9E77",
                  "#1B9E77","#1B9E77",
                  "#D95F02","#D95F02",
                  "#1B9E77","#1B9E77",
                  "#1B9E77","#1B9E77",
                  "#1B9E77","#1B9E77",
                  "#1B9E77","#1B9E77",
                  "#D95F02","#D95F02",
                  "#1B9E77","#1B9E77",
                  "#1B9E77","#1B9E77")

rv <- rowVars(assay(vsd))
# top 2000 variable genes 
select <- order(rv, decreasing=T)[seq_len(min(2000,length(rv)))] 
my_palette <- colorRampPalette(c("blue", "white", "red"))(1024)
heatmap3(assay(vsd)[select,], margins = c(10,5),
         hclustfun = hclust_comp, 
         distfun = dist_comp,
         col=my_palette,
         labRow = F,
         cexCol = 0.8,
         ColSideColors = colsidecolors,
         ColSideLabs = F, 
         balanceColor = F)

# 5 x 5.5

# saving 5 x 7
dev.copy(pdf, "2015-2-11-DESeq2-batch3-heatmap.pdf", width=5, height=5)
dev.off()

# select top 2000 gene nanes
heatmap_names <- rv[select]
heatmap_data <- merge(as.data.frame(heatmap_names), as.data.frame(resEvsM), by= 'row.names',sort=F)
rownames(heatmap_data) <- heatmap_data$Row.names
heatmap_data$Row.names <- NULL
heatmap_data$heatmap_names <- NULL
heatmap_data <- merge(heatmap_data, as.data.frame(counts(ddsEvsM,normalized=T)), by='row.names',sort=F)
write.csv(heatmap_data,"2015-2-11-heatmap3-2000genesvsd-Normalizedvalues.csv")

# PCA analysis
rv <- rowVars(assay(vsd))
sel <- order(rv, decreasing=T)[seq_len(min(2000,length(rv)))] # 500-2000 variable genes
pc <- prcomp(t(assay(vsd)[sel,]))

# 2D PCA plot with ggplot2
scores <- data.frame(sampleFiles, pc$x, sampleCons,  stringsAsFactors=T) # be sure to add grouping # SEE BELOW

my_palette <- c("#D95F02","#1B9E77")
# or factor(grouping)
(pcaplot <- ggplot(scores, aes(x=PC1, y=PC2, col=scores$sampleCons))) 
(pcaplot <- pcaplot + geom_point(size = 4) +
   ggtitle("Principal Components") +
   scale_color_manual(values = my_palette) + 
   theme(
  plot.title = element_text(face='bold'),
  legend.position="bottom", 
  legend.key = element_rect(fill='NA'),
  legend.title = element_blank(),
  legend.text = element_text(size=10, face="bold"),
  axis.text.y=element_text(colour="Black"), 
  axis.text.x=element_text(colour="Black"),
  axis.title.x=element_text(face="bold"), 
  axis.title.y=element_text(face='bold'),
  panel.grid.major.x=element_blank(), 
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x=element_blank(), 
  panel.grid.minor.y=element_blank(), 
  panel.background=element_rect(color='black',fill=NA)))
ggsave(pcaplot, file = "2015-2-11-batch3-2D-PCA_ggplot2_2group_bc.pdf", width = 4, height = 4)

# using grouping to ID sub classes
# set condition
condition <- factor(c("Ecad-","Ecad+", "Ecad-","Ecad+", 
               "Ecad-","Ecad+", "Ecad-","Ecad+", 
               "Ecad-","Ecad+", "Ecad-","Ecad+",
               "Ecad-","Ecad+", "Ecad-","Ecad+",
               "Ecad-","Ecad+", "Ecad-","Ecad+",
               "Ecad-","Ecad+", "Ecad-","Ecad+",))
condition <- as.integer(condition)

grouping <- factor(c("Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Non-Canonical","Non-Canonical"))


################################################################################
#
#
# Perform DE gene expression analysis between Canonical and Non-Canincal samples 
#
# 
#
################################################################################

ddsHTseqNCvsC <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                       directory = directory, 
                                       design= ~samplecons) 

colData(ddsHTseqNCvsC)$condition <- factor(colData(ddsHTseqNCvsC)$samplecons, 
                                      levels=c("Non-Canonical",'Canonical'))

# control for batch
ddsNCvsC <- DESeq(ddsHTseqNCvsC)
design(ddsNCvsC) <- formula(~ samplebatch + samplecons)
ddsNCvsC <- DESeq(ddsNCvsC)
resNCvsC <- results(ddsNCvsC)
resNCvsC <- resNCvsC [order(resNCvsC $padj),]
table(resNCvsC $padj < 0.1)
resdataNCvsC <- merge(as.data.frame(resNCvsC ), as.data.frame(counts(ddsNCvsC,normalized=T)), by='row.names',sort=F)
names(resdataNCvsC)[1] <- 'Gene'
head(resdataNCvsC)
write.csv(resdataNCvsC, file="2015-2-12-DESeq2-batch3-NCvsC-results-with-normalied.csv")

#subsetting on DE genes with padj < 0.1 and LFC > 1 and < -1.0
res.sub <- subset(resNCvsC, padj <  0.1) # subsetting DE genes with padj < 0.1 
res.sub <- subset(res.sub, log2FoldChange > 1.0 | log2FoldChange < -1.0) # subset on DE gene with FC greater than 1 
resdata <- merge(as.data.frame(res.sub), as.data.frame(assay(vsd)), by = 'row.names', sort=F)
names(resdata)[1] <- "gene"
write.csv(resdata, file = "2015-2-11-DESeq2-NCvsC-padj0.1-lfc-1-vsdMF-normalized.csv")

################################################################################
#
#
# Likelihood Ratio Test (LRT) in deseq2 to compare three group comparision 
#
# MAYBE NOT DO THIS FOR TIME BEING 
#
################################################################################

# setting condition (L1, L2, L3) for Lilkihood ratio test
ddsHTseqLRT <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                       directory = directory, 
                                       design=~condition)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseqLRT

# this may not owork
colData(ddsHTseqLRT)$condition<-factor(colData(ddsHTseqLRT)$condition, levels=c('L1','L2','L3'))

# gut of DESeq2 analysis
ddsLRT <- DESeq(ddsHTseqLRT, test='LRT', full=~cdh1_status+condition, reduced=~cdh1_status)

resultsNames(ddsLRT)

resLRT <- results(ddsLRT, name="condition_L2_vs_L1")
resLRT <- resLRT[order(resLRT$padj),]
head(resLRT)
resdat <- merge(as.data.frame(resLRT), as.data.frame(counts(ddsLRT,normalized=T)), by='row.names',sort=F)
names(resdat)[1] <- 'gene'
head(resdat)
write.csv(resdat, file="2015-2-11-results-LRT-DESeq2-with-normalized-batch3.csv")


################################################################################
#
#
# multifactor designs - controlling for three different batches of RNAseq 
# can analysis with more than one factor influencing the counts 
# from manual section 1.5
#
#
################################################################################

ddsMF <- dds
design(ddsMF) <- formula(~ Batch + condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
resMF <- resMF[order(resMF$padj),]
head(resMF)

write.csv(as.data.frame(resMF),file='2015-2-11-DESeq2-batch3_allsamples_withBATCHcontrol.csv')

table(resMF$padj < 0.1)
resdataMF <- merge(as.data.frame(resMF), as.data.frame(counts(ddsMF,normalized=T)), by='row.names',sort=F)
names(resdataMF)[1] <- 'gene'
head(resdataMF)
write.csv(resdataMF, file="2015-2-11-DESeq2-batch3-final-results-with-normalied.csv")

# variance stabilization normalization
vsdMF <- varianceStabilizingTransformation(ddsMF,blind=T)
# clustering analysis
distsRL <- dist(t(assay(vsdMF)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2015-2-11-deseq2_heatmaps_withoutPD8640_PD1849.png")
dev.off()

plotMA(ddsMF, ylim=c(-8,8),main = "Ecad+/- RNAseq batch affect controlled")
dev.copy(png, "2015-2-11-DESeq2_MAplot_batchcontrolled_allsamples2.png")
dev.off()



# views clustering on individual datasets in unbiased way
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(ddsNCvsC)[,1] / sizeFactors(ddsNCvsC)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"2015-2-11-DESeq2_variance_stabilizing.png")
# axis is square root of variance over the mean for all samples
dev.off()

# model-based clustering alorithm (mclust)
# selecting model - multivariat mixture = VII = diagonal, varying volume and shape
fit <- Mclust(pc$x, G=3, modelNames = "VII")
fit$classification
# 3D pca with rgl plot3d using same pc 
plot3d(pc$x,size=15,lit=T, box=F, alpha=0.75, axes=F, col=grouping)
text3d(pc$x, text=sampleFiles, adj=-0.2, font=1, cex=0.8, main="Principal Components")
grid3d(side='z',at=list(z=0))
# plot3d(pc$x, col=pca$cluster)
# save as postscript pdf file then edit in gimp
rgl.snapshot("2015-2-11-3D-pca_view2.png",'png')

# or use the pca3d package 
# I prefere this 
library("pca3d")
threedpca <- pca3d(pc, show.plane=F, col=scores$grouping, radius=2, axes.color='black',
      show.axe.titles=F)
rgl.snapshot("2014-10-10-pca3d.grouping.png",'png')
save(threedpca, file="2014-10-10-3d-PCA_3group.RData")

# optional K-means based clustering
km <- kmeans(assay(ddsNCvsC),3) # not sure how to integrate into plot
pc$cluster <- as.integer(km$cluster)
plot3d(pc$x,size=10,lit=T, axes=F,box=F)
text3d(pc$x, text=sampleFiles, adj=-0.2, font=0.5, col=colorRampPalette(c("red",'red'))(2))
plot3d(pc$x, col=pca$cluster)

save.image("~/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_batch2/2014-10-10-batch2.RData")

#### replacing outlier value with estimated value as predicted
#### by distrubution using "trimmed mean" approach
#### recommended if you have several replicates per treatment
#### DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file='2014-6-10-replaceoutliers-results1.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file='2014-6-10-replaceoutliers-results2.csv')
