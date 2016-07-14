\name{combineDesignMatrix}
\alias{combineDesignMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Combining RNA- and Ribo-seq treatment conditions into one single
  design matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combineDesignMatrix(rnaCond, riboCond)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rnaCond}{
%%     ~~Describe \code{rnaCond} here~~
}
  \item{riboCond}{
%%     ~~Describe \code{riboCond} here~~
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
function (rnaCond, riboCond) 
{
    if (!is.data.frame(rnaCond)) 
        rnaCond <- data.frame(cond = rnaCond)
    if (!is.data.frame(riboCond)) 
        riboCond <- data.frame(cond = riboCond)
    if (ncol(rnaCond) != ncol(riboCond)) 
        stop("rna- and ribo-seq data must have the same number of conditions")
    numCond <- ncol(rnaCond)
    numRNASmps <- nrow(rnaCond)
    numRiboSmps <- nrow(riboCond)
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
    extendedConds <- paste0("combinedCond$", colnames(combinedCond))
    fmla <- as.formula(paste("~", paste(extendedConds, collapse = "+")))
    model.matrix(fmla)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
