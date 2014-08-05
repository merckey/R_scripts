# DESeq2 run Benjamini & Hockberg (BH) to control the false discovery rate
# use reported p values to calculate fdr rate (which is what NYU did) 
# fdr actually lowered the number of significant DE genes

# produce data frame 'res' using DESeq2 pipeline
res <- read.csv(file='EcadminusVSplusgene_exp-diff.csv')

head(res)

# should have pvalue and padj fields (this is BH)
# to calculate fdr instead

p <- res[,12]
# takes pvalues from 5th colume of data frame 'res'
length(p)

# ?p.adjust
# if 'none' is selected for methods = it just prints unadjusted pvals

fdr <- p.adjust(p,'fdr', length(p))

head(fdr)


# combine newly calculated fdr rates into a new Res data frame
# this crashed the R session last time.
res.p <- data.frame(res, fdr)
res.p <- res.p[order(res.p$fdr),]

# write file to csv!
write.csv(as.data.frame(res.p),file='corrected_NYU_fdr_results.csv')
