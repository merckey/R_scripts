# DESeq2 run Benjamini & Hockberg (BH) to control the false discovery rate
# use reported p values to calculate fdr rate (which is what NYU did) 
# fdr is al ittle more flexibile in finding significance (not as stringent)

# produce data frame 'res' using DESeq2 pipeline

head(res)

# should have pvalue and padj fields (this is BH)
# to calculate fdr instead

p <- res[,5]
# takes pvalues from 5th colume of data frame 'res'

fdr <- p.adjust(p, p.adjust.methods, n = length(p))

head(fdr)


# combine newly calculated fdr rates into a new Res data frame
# this crashed the R session last time.
res.p <- data.frame(res, fdr)
res.p <- res.p[order(res.p$fdr),]

# write file to csv!
write.csv(as.data.frame(res.p),file='corrected_fdr_results.csv')
