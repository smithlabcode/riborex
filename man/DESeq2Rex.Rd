\name{DESeq2Rex}
\alias{DESeq2Rex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DESeq2Rex(rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, minMeanCount = 1)
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
        stop("rna- and ribo-seq data must have the same set of genes")
    if (ncol(rnaCond) != ncol(riboCond)) 
        stop("rna- and ribo-seq data must have the same number of conditions")
    if (!is.data.frame(rnaCond)) 
        rnaCond <- data.frame(cond = rnaCond)
    if (!is.data.frame(riboCond)) 
        riboCond <- data.frame(cond = riboCond)
    keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
    keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
    keep <- intersect(keep.rna, keep.ribo)
    rnaCntTable <- rnaCntTable[keep, ]
    riboCntTable <- riboCntTable[keep, ]
    numCond <- ncol(rnaCond)
    numRNASmps <- nrow(rnaCond)
    numRiboSmps <- nrow(riboCond)
    combCntTbl <- cbind(rnaCntTable, riboCntTable)
    rnaExpansion <- matrix(factor(rep(rep(0, numCond + 1), numRNASmps)), 
        nrow = numRNASmps)
    rnaCond <- cbind(rnaCond, as.data.frame(rnaExpansion))
    numExtendedCols <- length(colnames(rnaCond))
    colnames(rnaCond)[(numCond + 1):numExtendedCols] <- paste0("extra", 
        seq(numCond + 1))
    riboCond <- cbind(riboCond, intercept = factor(1), riboCond)
    colnames(riboCond)[(numCond + 1):numExtendedCols] <- paste0("extra", 
        seq(numCond + 1))
    combinedCond <- rbind(rnaCond, riboCond)
    extendedConds <- colnames(combinedCond)
    fmla <- as.formula(paste("~", paste(extendedConds, collapse = "+")))
    dds <- DESeqDataSetFromMatrix(countData = combCntTbl, colData = combinedCond, 
        design = fmla)
    dds <- DESeq(dds)
    if (is.null(contrast)) {
        res <- results(dds)
    }
    else {
        res <- results(dds, contrast = contrast)
    }
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
