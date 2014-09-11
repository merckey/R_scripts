################################
# 2014-8-13
# Comparrison of EMT-noncanonical samples with PD 9210 includes
# will do calculations for batch control

library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_Canonical_vs_Non-Canonical/non-Canonical/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_Canonical_vs_Non-Canonical/non-Canonical/'

sampleFiles <- grep('counts',list.files(directory),value=T)

# view sampleFiles
sampleFiles


sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch1","Batch1","Batch2","Batch2")

# set sampleConditions for Nicole's RNAseq ecad neg/pos experiment
sampleCondition <- c('E_minus','E_plus')
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition,
                          Batch = sampleBatch)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory = directory, 
                                       design=~condition)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, 
                                    levels=c('E_plus','E_minus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)

# multifactor designs
# can analysis with more than one factor influencing the counts 
# from manual section 1.5

ddsMF <- dds

design(ddsMF) <- formula(~ Batch + condition)

ddsMF <- DESeq(ddsMF)

plotMA(ddsMF, ylim=c(-8,8),main = "Ecad+/- RNAseq batch effect controlled")
dev.copy(png, "2014-8-12-DESeq2_MAplot_batchcontrolled_final.png")
dev.off()


# did not do this 2014-8-13
rld <- rlogTransformation(ddsMF, blind=T)
#print(plotPCAWithSampleNames(rld2, intgroup=c('condition')))
#dev.copy(png, "2014-8-12-DESeq2_MAplot_PCA_batchcontrol_final.png")
#dev.off()

resMF <- results(ddsMF)

resMF <- resMF[order(resMF$padj),]

head(resMF)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# padj should be ranked lowest adj pval to high (most sig to least sig)
# save data 'res' to csv!
write.csv(as.data.frame(resMF),file='2014-8-13-DESeq2_nonCanonical_Batchcontrol_final.csv')

# clustering analysis
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsMF),paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-8-12-DESeq2_batchcontrol_clustering_final.png")
dev.off()
