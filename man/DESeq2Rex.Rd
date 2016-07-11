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
DESeq2Rex(rnaCntTable, riboCntTable, rnaCondition, riboCondition)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rnaCntTable}{
%%     ~~Describe \code{rnaCntTable} here~~
}
  \item{riboCntTable}{
%%     ~~Describe \code{riboCntTable} here~~
}
  \item{rnaCondition}{
%%     ~~Describe \code{rnaCondition} here~~
}
  \item{riboCondition}{
%%     ~~Describe \code{riboCondition} here~~
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
function (rnaCntTable, riboCntTable, rnaCondition, riboCondition) 
{
    library(DESeq2)
    combCntTbl <- cbind(rnaCntTable, riboCntTable)
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps), 
        rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(RNACntTable)), rep(1, ncol(RiboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(RNACntTable) + numCtlRiboSmps), 
        rep(1, numCaseRiboSmps)))
    combColData <- data.frame(condition = condition)
    combColData$dataType <- dataType
    combColData$RiboDiff <- RiboDiff
    rownames(combColData) <- colnames(combCntTbl)
    dds <- DESeqDataSetFromMatrix(countData = combCntTbl, colData = combColData, 
        design = ~condition + dataType + RiboDiff)
    dds <- DESeq(dds)
    res <- results(dds)
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
