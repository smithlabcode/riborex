\name{edgeRDRex}
\alias{edgeRDRex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edgeRDRex(rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, minMeanCount = 1)
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
    dge.rna <- DGEList(counts = rnaCntTable)
    dge.rna <- calcNormFactors(dge.rna)
    design.rna <- dataFrameToDesignMatrix(rnaCond)
    dge.rna <- estimateDisp(dge.rna, design.rna)
    dge.ribo <- DGEList(counts = riboCntTable)
    dge.ribo <- calcNormFactors(dge.ribo)
    design.ribo <- dataFrameToDesignMatrix(riboCond)
    dge.ribo <- estimateDisp(dge.ribo, design.ribo)
    dispersion.rna <- getDispersion(dge.rna)
    dispersion.ribo <- getDispersion(dge.ribo)
    dispersion <- matrix(c(rep(dispersion.rna, ncol(rnaCntTable)), 
        rep(dispersion.ribo, ncol(riboCntTable))), nrow = dim(rnaCntTable)[1])
    combCntTbl <- cbind(rnaCntTable, riboCntTable)
    combFactors <- c(dge.rna$samples$norm.factors, dge.ribo$samples$norm.factors)
    dge <- DGEList(counts = combCntTbl, norm.factors = combFactors)
    design <- combineDesignMatrix(rnaCond, riboCond)
    fit <- glmFit(dge, design = design, dispersion = dispersion)
    lrt <- glmLRT(fit, contrast = contrast)
    topGenes <- topTags(lrt, n = Inf)
    topGenes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
