\name{edgeRRex}
\alias{edgeRRex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edgeRRex(rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, minMeanCount = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rnaCntTable}{
%%     ~~Describe \code{rnaCntTable} here~~
}
  \item{riboCntTable}{
%%     ~~Describe \code{riboCntTable} here~~
}
  \item{rnaCond}{
%%     ~~Describe \code{rnaCond} here~~
}
  \item{riboCond}{
%%     ~~Describe \code{riboCond} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{minMeanCount}{
%%     ~~Describe \code{minMeanCount} here~~
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
function (rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, 
    minMeanCount = 1) 
{
    if (!identical(rownames(rnaCntTable), rownames(riboCntTable))) 
        stop("RNA- and Ribo-seq data must have the same set of genes")
    if (ncol(rnaCond) != ncol(riboCond)) 
        stop("RNA- and Ribo-seq data must have the same number of conditions")
    if (!is.data.frame(rnaCond)) 
        rnaCond <- data.frame(cond = rnaCond)
    if (!is.data.frame(riboCond)) 
        riboCond <- data.frame(cond = riboCond)
    keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
    keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
    keep <- intersect(keep.rna, keep.ribo)
    rnaCntTable <- rnaCntTable[keep, ]
    riboCntTable <- riboCntTable[keep, ]
    combCntTbl <- cbind(rnaCntTable, riboCntTable)
    dge <- DGEList(counts = combCntTbl)
    dge <- calcNormFactors(dge)
    design <- combineDesignMatrix(rnaCond, riboCond)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, contrast = contrast)
    topGenes <- topTags(lrt, n = Inf)
    topGenes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
