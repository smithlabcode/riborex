\name{voomRex}
\alias{voomRex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Identifying differential translation from ribosome profiling data 
  using Voom as an engine
}
\description{
  This function identifies differentially translated genes by
  co-analyzing both RNA- and Ribo-seq data. It uses Voom
  as its engine.
}
\usage{
voomRex(rnaCntTable, riboCntTable, rnaCond, riboCond, contrast = NULL, minMeanCount = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rnaCntTable}{
  a matrix or data frame of raw read counts summarized from RNA-seq data.
  The rows correspond to genes and columns correspond to samples.
}
  \item{riboCntTable}{
  a matrix or data frame of raw read counts summarized from Ribo-seq data.
  The rows correspond to genes and columns correspond to samples.
  It should have the same set of genes as the rnaCntTable.
}
  \item{rnaCond}{
  a vector for single factor experiment or a data frame for multi-factor
  experiment. The order of samples should be the same as in rnaCntTable.
}
  \item{riboCond}{
  a vector for single factor experiment or a data frame for multi-factor
  experiment. The order of samples should be the same as in riboCntTable.
  The set of factors and the order of them should be the same as rnaCond.
}
  \item{contrast}{
  a vector to specify which contrast to be applied on.
}
  \item{minMeanCount}{
  The criterion for filtering low read count. The genes whose mean read
  count across all RNA- and Ribo-seq samples less than minMeanCount
  will be omitted.
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
  Wenzheng Li <wenzhenl@usc.edu>
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
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    if (!is.null(contrast)) {
        fit <- contrasts.fit(fit, contrasts = contrast)
    }
    fit <- eBayes(fit)
    topGenes <- topTable(fit, coef = ncol(design), number = Inf)
    topGenes
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
