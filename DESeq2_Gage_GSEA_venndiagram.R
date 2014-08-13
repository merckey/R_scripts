# 

library(limmma)

canonical <- read.delim("~/Desktop/RNAseq_Nicole_Ecad/Gage_pathview/Canonical_group/2014-8-6-Canoical-KEGG.txt",sep='\t')
noncanonical <- read.delim("~/Desktop/RNAseq_Nicole_Ecad/Gage_pathview/Non-canonical_group/2014-8-6-Noncanoical-KEGG.txt",sep='\t')

can.sig <- canonical$p.geomean[1:31]
noncan.sig <- noncanonical$p.geomean[1:33]

merged <- cbind(can.sig,noncan.sig)

counts <- matrix(0, nrow=length(merged),ncol=2)
counts <- vennCounts(merged)

counts
