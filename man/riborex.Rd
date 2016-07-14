\name{riborex}
\alias{riborex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Identifying differential translation from ribosome profiling data 
}
\description{
  This function identifies differentially translated genes by
  co-analyzing both RNA- and Ribo-seq data. It uses existing
  RNA-seq data analysis tools as its engines. Currently, the
  engines include DESeq2, edgeR and Voom. With edgeR, you can
  choose estimate dispersion of RNA- and Ribo-seq data separately
  or together.
}
\usage{
riborex(rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, minMeanCount = 1, engine = "DESeq2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rnaCntTable}{
  A data.frame.
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
  \item{engine}{
%%     ~~Describe \code{engine} here~~
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
    minMeanCount = 1, engine = "DESeq2") 
{
    if (!identical(rownames(rnaCntTable), rownames(riboCntTable))) 
        stop("RNA- and Ribo-seq data must have the same set of genes")
    if (ncol(rnaCond) != ncol(riboCond)) 
        stop("RNA- and ribo-seq data must have the same number of conditions")
    if (engine == "DESeq2") {
        DESeq2Rex(rnaCntTable, riboCntTable, rnaCond, riboCond, 
            contrast, minMeanCount)
    }
    else if (engine == "edgeR") {
        edgeRRex(rnaCntTable, riboCntTable, rnaCond, riboCond, 
            contrast, minMeanCount)
    }
    else if (engine == "edgeRD") {
        edgeRDRex(rnaCntTable, riboCntTable, rnaCond, riboCond, 
            contrast, minMeanCount)
    }
    else if (engine == "Voom") {
        voomRex(rnaCntTable, riboCntTable, rnaCond, riboCond, 
            contrast, minMeanCount)
    }
    else {
        stop("Error: unrecognized engine name")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
