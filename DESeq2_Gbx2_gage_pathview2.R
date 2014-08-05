# gage pathway analysis of GBX2-high group with Intersection not empty specification for summarize overlaps



setwd("~/Desktop/RNAseq_Nicole_Ecad/Gage_pathview/")
library("BiocInstaller", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
biocLite(c('pathview','gage','gageData','GenomicAlignments','TxDb.Mmusculus.UCSC.mm10.knownGene'))
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
exByGn <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, 'gene')

library(Rsamtools)
fls <- list.files(pattern='bam$',full.names=T)
bamfls <- BamFileList(fls)
flag <- scanBamFlag(isNotPrimaryRead=F, isProperPair=T)
param <- ScanBamParam(flag=flag)
gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="IntersectionNotEmpty",ignore.strand=T,single.end=F, param=param)
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

library(DESeq2)
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
fc.kegg.p <- gage(exp.fc, gsets = kegg.sigmet, ref =NULL, samp=NULL)
sel <- fc.kegg.p$greater[,'q.val'] < 0.1 &
  !is.na(fc.kegg.p$less[,'q.val'])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  +            !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:10],
                      function(pid) pathview(gene.data=exp.fc, 
                                             pathway.id = pid,
                                             species = "mouse", 
                                             out.suffix=out.suffix))

# to see axon guidance genes
# > kegg.g2$kg.sets$mmu04360
