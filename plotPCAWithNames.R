plotPCAWithSampleNames = function(x, intgroup="condition", ntop=200)  # can change ntop from 200-500
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