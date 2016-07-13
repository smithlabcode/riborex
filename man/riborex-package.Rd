\name{riborex-package}
\alias{riborex-package}
\alias{riborex}
\docType{package}
\title{
    A R package for identification of differential translation from Ribo-seq data
}
\description{
\packageDescription{
    Riborex provides functions to identify differentially translated genes
    from Ribo-seq data. Riborex is built upon existing RNA-seq data analysis
    tool. Currently, we implemented riborex based on DESeq2, edgeR and Voom.
    For edgeR, a version in which dispersion of RNA- and Ribo-seq data
    are estimated independently is also implemented.
  }
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{riborex}
\packageIndices{riborex}
functions:

riborex
DESeq2Rex
edgeRRex
edgeRDRex
voomRex
}
\author{
\packageAuthor{
  Wenzheng Li, Weili Wang and Andrew D. Smith
}

Maintainer: \packageMaintainer{Wenzheng Li <wenzhenl@usc.edu>}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
