\name{plotPCAWithSampleNames}
\alias{plotPCAWithSampleNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotPCAWithSampleNames(x, intgroup = "condition", ntop = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{intgroup}{
%%     ~~Describe \code{intgroup} here~~
}
  \item{ntop}{
%%     ~~Describe \code{ntop} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, intgroup = "condition", ntop = 500) 
{
    library("genefilter")
    library("lattice")
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    names = colnames(x)
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]), 
        1, paste, collapse = " : "))
    if (nlevels(fac) >= 3) 
        colours = brewer.pal(nlevels(fac), "Dark2")
    else colours = c("lightgreen", "dodgerblue")
    xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
        pch = 16, cex = 1.5, panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            ltext(x = x, y = y, labels = names, pos = 1, offset = 0.8, 
                cex = 0.7)
        }, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
            text = list(levels(fac)), rep = FALSE)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
