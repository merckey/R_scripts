### bioinformatic analysis of RNA-seq data
### recipe is from Dave Wheeler's blog at Massey University
### http://dwheelerau.com/2013/04/15/how-to-use-deseq-to-analyse-rnaseq-data/
###
###
### Should have merged counts file from HT-seq-count
### be sure to replace "experimentname" with actual name of experiment 
### be sure to set working directory to where data files are located

library("DESeq")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_GBX2_group/")
countsTable <- read.delim('merged_counts.txt',header=T,row.names='gene')
cou#rownames(countsTable) <- countsTable$gene 
countsTable <- countsTable[,1]

### add experimental treatment/condition information for each column
# per example
conds  <- factor(c("Eplus","Eplus","Eplus",'Eminus','Eminus','Eminus'))

cds <- newCountDataSet(data, conds)

cds <- estimateSizeFactors(cds)

### to visualize size favtors
sizeFactors(cds)

## differental expression calculation 
cds <- estimateDispersions(cds)

## plot dispersion estimate
plotDispEsts <- function(cds){
  plot(rowMeans(counts(cds,normalized=T)), fitInfo(cds)
  $perGeneDispEsts, pch = '.',log='xy',ylab='dispersion',
  xlab='mean of normalized counts')
  xg = 10^seq(-.1,5,length.out=300)
  lines(xg, fitInfo(cds)$dispFun(xg),col='red')
}
plotDispEsts(cds)
jpeg("2014-7-23-DispEsts_DESeq.jpg")
plotDispEsts(cds)

##################################################################
# I am adding Dave Wheeler's explaination of dispersion here 
#
# Remember that the level of dispersion is related to the 
# biological variation seen in each treatment. Note this is a 
# double log graph, ahhh my brain hurts, but look carefully at 
# what is happening on the Y axis, as read count increases 
# dispersion decreases, as we would expect. However, a lot of 
# those black dot values below the red fitted line are probably 
# underestimates of dispersion based on the small samples sizes 
# used (only 3 replicates values measured in this example). 
# So, to be conservative DESeq moves all these values below the 
# red line up to the fitted value, BUT keeps all the empirical 
# values above the red line, even if some of these are over estimates 
# of dispersion, as I said it is being conservative. Why does 
# dispersion matter? The more dispersion, or biological variation, 
# the bigger the difference between counts between treatments is 
# required before it differences become significant!
###################################################################
       
# add differential expression values
# make sure to change sample names to actual experimental conditions
res <- nbinomTest(cds, "Eplus",'Eminus')

# visualize res
head(res)

jpeg("plotMA.jpg")
plotMA(res)
def.off

# compare normalized counts between replicates 
## besure to change name of "untreated" to actual experimetnal sample
ncu <- counts(cds, normalized=T)[, conditions(cds)=='untreated']

jpeg("MA_untreated_only.jpg")
plotMA(data.frame(baseMean = rowMeans(ncu),log2FoldChange = log2( ncu[,2] / ncu[,1] )),
       col = 'black')
def.off()

# make histogram of p value distrubtion
jpeg("hist_pval.jpg")
hist(res$pval, breaks=100, col='skyblue',border='slateblue',main='')
dev.off()

# create subset table base on FDR cut offs correcting for multiple hypotheses
# setting FDR at 10% for now
res <- res[order(res$padj),]
# visualize FDR
head(res)

# save data!
write.csv(res, file="Experimentalname.csv")

# save FDR data!
write.csv(resSig, file='Experimentalname_0.01FDR.csv')


##################################################################
# I am again, adding Dave Wheeler's explaination of process so far
# this *should* conclude basic analysis
# 
# We are done with the basic analysis, we can do some more 
# sophisticated filtering of the data. The main goal is to remove 
# tests that give little chance of producing a meaningful result, 
# with the aim of keeping the FDR in check. This results in an 
# increased power to detect significant results whilst keeping 
# the same level of type I error the same (false positive rate)
###################################################################

# take sum of counts across all exp conditions
rs <- rowSums(counts(cds))

# filter out the lowest 40% quantile
theta <- 0.4

use <- (rs > quantile(rs, probs = theta))

table(use)

cdsFilt <- cds[ use, ]

# run binomal test on filtered data
resFilt <- nbinomTest(cdsFilt, "Experimentalcondition1", 'Experimentalcondition2')

# look at new FDR
tabAll <- table(res$padj<0.1)
tabFilt <- table(resFilt$padj<0.1)

tabAll
tabFilt

# you are basically removing NA values in above filtering

# plot barplot of filterd pvals
# powderblue are pvals that passed filtering
h1 <- hist(res$pval,breaks=50,plot=F)
h2 <- hist(resFilt$pval, breaks=50,plot=F)
colori <- c("do not pass"='khaki','pass'='powderblue')
barplot(height=rbind(h1$counts,h2$counts),beside=F,col=colori,
        space=0, main ='',ylab='frequency')
text(x=c(0,length(h1$counts)),y=0,label=paste(c(0,1)),adj=c(0.5,1.7), xpd=NA)
legend("topright",fill=rev(colori),legend=rev(names(colori)))
dev.copy(png,'hist_filthist.png')
def.off()


cdsBlind <- estimateDispersions(cds, method='blind')

vsd <- varianceStabilizingTransformation(cdsBlind)

# plot per gene stdev against rank of the mean, stdev of transfromated data
# across samples against hte mean 
souce("http://bioconductor.org/biocLite")
biocLite("vsn")
library('vsn')

par(mfrow=c(1,2))

notAllZero <- (rowSums(counts(cds))>0)

meanSdPlot(log2(counts(cds)[notAllZero, ] +1), ylim = c(0,2.5))

meanSdPlot(vsd[notAllZero, ], ylim = c(0,2.5))

mod_lfc <- (rowMeans(exprs(vsd)[ , conditions(cds)=='experimentalconditionname1',drop=F]) -
              rowMeans(exprs(vsd)[ , conditions(cds)=='experimentalconditionname2', drop=F]))

lfs <- res$log2FoldChange

table(lfcs[!is.finite(lfc)],useNA='always')

logdecade <- 1 + round(log10(1+rowMeans(counts(cdsBlind,normalized=T))))

lfccol <- colorRampePalette(c('gray','blue'))(6)[logdecade]

ymax=5

plot(pmax(-ymax,pmin(ymax,lfc)),mod_lfc,
     xlab = 'ordinary log-ratio', ylab = 'moderated log-ratio',
     cex = 0.45, asp = 1, col = lfccol,
     pch = ifelse(lfcs<(-ymax), 60, ifelse(lfc>ymax, 62, 16)))

abline(a=0, b=1, col='red3')

# graph should show weakly to strongly exp genes in grey to blue

# produce heat maps to see patterns using transformed log data
library("gplots")
library('RColorBrewer')

cdsBlind <- estimateDispersions(cds, method='blind')
vsd <- varianceStabilizingTransformation(cdsBlind)
select <- order(rowMeans(counts(cds)), decreasing =F)[1:250]
colfunc <- colorRampPalette(c("blue",'red'))
heatmap.2(exprs(vsd)[select,],col=colfunc(100), trace = 'none', margin = c(10,6))

# run principal component analysis on data
# good for visualizing effect of experimental covariats and batch effect
# ideal for examining primary and matching mets
print(plotPCA(vsd, intgroup=c("condition")))

sessionInfo()
