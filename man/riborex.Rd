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
  engines include DESeq2, edgeR, edgeRD and Voom.
  Among them, DESeq2, edgeR and Voom estimate dispersion across
  all RNA- and Ribo-seq samples, while edgeRD estimates dispersion
  for RNA- and Ribo-seq data separately.
}
\usage{
riborex(rnaCntTable, riboCntTable, rnaCond, riboCond, engine = "DESeq2", contrast = NULL, minMeanCount = 1)
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
  experiment to label the treatment applied on samples. The order of
  samples should be the same as that in rnaCntTable.
}
  \item{riboCond}{
  a vector for single factor experiment or a data frame for multi-factor
  experiment to label the treatment applied on samples. The order of
  samples should be the same as that in riboCntTable.
  The set of factors and the order of them should be the same as rnaCond.
}
  \item{engine}{
  A string to specify the name of engine to be used. You can choose
  "DESeq2", "edgeR", "edgeRD" or "Voom". Among them, DESeq2, edgeR
  and Voom estimate dispersion across
  all RNA- and Ribo-seq samples, while edgeRD estimates dispersion
  for RNA- and Ribo-seq data separately.
}

  \item{contrast}{
    a vector to specify the comparison to extract the results.
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

\examples{
  ### load the data
  data(riborexdata)
  ### get rna-seq read count table
  rna <- riborexdata$rna
  ### get ribo-seq read count table
  ribo <- riborexdata$ribo
  ### prepare rna-seq condtions
  rnacond <- factor(c("control", "control", "treated", "treated"))
  ### prepare ribo-seq condtions
  ribocond <- factor(c("control", "control", "treated", "treated"))
  ### run riborex with default engine "DESeq2"
  res.deseq2 <- riborex(rna, ribo, rnacond, ribocond)
  ### run riborex with engine "edgeR"
  res.edgeR <- riborex(rna, ribo, rnacond, ribocond, "edgeR")
  ### run riborex with engine "edgeRD"
  res.edgeRD <- riborex(rna, ribo, rnacond, ribocond, "edgeRD")
  ### run riborex with engine "Voom"
  res.Voom <- riborex(rna, ribo, rnacond, ribocond, "Voom")
}
