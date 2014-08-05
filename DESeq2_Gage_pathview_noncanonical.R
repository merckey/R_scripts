# gage pathway analyssi for non-canonical EMT samples:
# PD2329, PD2342 and PD798 - not using PD 2523 for initial pathview analysis
# as of 2014-8-6

setwd("~/Desktop/RNAseq_Nicole_Ecad/Gage_pathview/Non-canonical_group/")
#library("BiocInstaller", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
#biocLite(c('pathview','gage','gageData','GenomicAlignments','TxDb.Mmusculus.UCSC.mm10.knownGene'))
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(gage)
library(DESeq2)
library(Rsamtools)
library(org.Mm.eg.db)
library(pathview)
library(DEXSeq)
library(GenomicAlignments) # I think I only need this package and not 'Rsamtools'

MulticoreParam(workers=12)

exByGn <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, 'gene')

fls <- list.files(pattern='bam$',full.names=T)

# check Bam files are correct 
fls 

bamfls <- BamFileList(fls)
flag <- scanBamFlag(isNotPrimaryRead=F, isProperPair=T)
param <- ScanBamParam(flag=flag)
gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="IntersectionNotEmpty",
                           ignore.strand=T,single.end=F, 
                           param=param)
cnts=assay(gnCnt)
dim(cnts)

sel.rn = rowSums(cnts) != 0 
cnts = cnts[sel.rn,]
dim(cnts)
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
cnts.norm = log2(cnts.norm+8)
range(cnts.norm)

# DESeq2 analysis
grp.idx <- rep(c('Ecad_plus','Ecad_minus','Ecad_plus','Ecad_minus','Ecad_plus','Ecad_minus'))
coldata = DataFrame(grp=factor(grp.idx))
coldata
dds <- DESeqDataSetFromMatrix(cnts, colData=coldata, design = ~grp)
dds <- DESeq(dds)
deseq2.res <- results(dds)
deseq2.res2 <- deseq2.res[order(deseq2.res$padj),]
head(deseq2.res2)
deseq2.fc <- deseq2.res$log2FoldChange
names(deseq2.fc) = rownames(deseq2.res)
exp.fc = deseq2.fc 
out.suffix = 'deseq2'

require(gage)
kegg.g2 <- kegg.gsets(species = 'mmu',id.type='kegg')
kegg.sigmet <- kegg.g2$kg.sets[kegg.g2$sigmet.idx]
fc.kegg.p <- gage(exp.fc, gsets = kegg.sigmet, 
                  ref =NULL, samp=NULL)
head(fc.kegg.p$greater[,1:5],10)

# or use z-test 
fc.kegg.2p <- gage(exp.fc, gsets = kegg.sigmet, 
                   ref=NULL, samp=NULL, saaTest = gs.zTest)
head(fc.kegg.2p$greater[,1:5],10)

# 14-8-6 using z-test statistics - may have more false positives 

# can alter q.val from 0.1 for pathways passing FDR or just download everything by setting to 0.9
sel <- fc.kegg.2p$greater[,'q.val'] < 1.0 & !is.na(fc.kegg.p$less[,'q.val'])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.2p$less[, "q.val"] < 1.0 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.2p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

require(pathview)
#view first 10 pathways as demo
pv.out.list <- sapply(path.ids2[1:10],
                      function(pid) pathview(gene.data=exp.fc, 
                                             pathway.id = pid,
                                             species = "mouse", 
                                             out.suffix=out.suffix))

# graphviz 
pv.out.list2 <- sapply(path.ids2[1:10],
                      function(pid) pathview(gene.data=exp.fc, 
                                             pathway.id = pid,
                                             species = "mouse", 
                                             out.suffix=out.suffix,
                                             kegg.native=F,
                                             sign.pos='bottomleft'))

write.table(fc.kegg.2p$greater,file='2014-8-6-Noncanonical-greater.txt')
write.table(fc.kegg.2p$greater,file='2014-8-6-Noncanonical-greater.txt',sep='\t')
write.table(rbind(fc.kegg.2p$greater,fc.kegg.2p$lesser),file='2014-8-6-Noncanoical-KEGG.txt',sep='\t')