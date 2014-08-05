#####
# DEXSeq analysis for differential exon usage in RNA-seq data
# all samples - no canonical versus non-canonical

library("DEXSeq")
setwd("/Users/dballi/Desktop/RNAseq_Nicole_Ecad/DEXSeq_analysis")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/DEXSeq_analysis'

sampleFiles <- grep('sj.txt',list.files(directory), value=T)

sampleFiles
# should be list of each splice junction txt file e.g PD2204Eplus_sj.txt

flattenedFile <- grep('gff',list.files(directory), value=T)

flattenedFile
# should be the name of the gff anotation file 

# creat sample table of sample - condition - library type 
# library type for us isn't important but you can analyze single and paired end reads at same time
sampleTable <- data.frame( 
  row.names = c("PD1849Eplus","PD2204Eplus",'PD2329Eplus',"PD2342Eplus","PD2523Eplus","PD798Eplus",
                "PD1849Eminus","PD2204Eminus","PD2329Eminus",'PD2342Eminus','PD2523Eminus','PD798Eminus'),
  condition = c("Ecad_plus","Ecad_plus","Ecad_plus","Ecad_plus","Ecad_plus","Ecad_plus",
                "Ecad_minus","Ecad_minus","Ecad_minus","Ecad_minus","Ecad_minus","Ecad_minus"),
  libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end",
              "paired-end","paired-end","paired-end","paired-end","paired-end","paired-end"))
  
sampleTable  


## create DEXSeqDataSet
dxd = DEXSeqDataSetFromHTSeq(
  sampleFiles,
  sampleData=sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile
  )

# check dxd data
colData(dxd)
head(featureCounts(dxd),5)
head(rowData(dxd), 3 )

# normalisation 
dxd <- estimateSizeFactors( dxd )

# set processor cores to 4
BAPPARAM <- MulticoreParam(workers=4)

# perform standard differential exon usage analysis 
# this is taking a long time 
# be sure to set processor cores to more than 1
dxr <- DEXseq( dxd )

# to see description of layout of dxr1
elementMetadata(dxr)$description

# how many genes are FDR < 0.1
table(dxr$padj, 0.1)

# how many genes are affected
table(tapply(dxr$padj, 0.1, dxr$groupID,any))

# make MA plot - sig genes should be red in color
plotMA(dxr, cex=0.4)
dev.copy(png,'2014-7-1-DEXSeq-MAplot')
dev.off()

# order by FDR for most differentially expressed exons 
res <- dxr[order(dxr$padj),]
head(res)
write.csv(as.data.frame(res),file='2014-7-1-DEXSeq-analysis-pooledsamples.csv')

# order by log2fold_Ecad_minus_Ecad_plus
# there is a lot of func popping up bc somethings are just not expressed
# in one group versus the other 
# I would suggest stopping with sorting by padj
res2 <- dxr[order(dxr$log2fold_Ecad_minus_Ecad_plus),]
head(res2)
write.csv(as.data.frame(res2),file='2014-7-1-DEXSeq-analysis-pooledsamples_ordered_by_log2.csv')

# plot fitted expression for Ncam1
plotDEXSeq(dxr, "chr3_S100a7a+",legend=T,cex.axis=1.2,cex=1.3,lwd=2)
dev.copy(png, '2014-7-7-DEXSeq-s100a7a')
dev.off()

