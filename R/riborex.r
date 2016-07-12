# Copyright (C) 2016 University of Southern California and
#                    Wenzheng Li, Weili Wang and  Andrew D. Smith
#
# Authors: Wenzheng Li and Weili Wang
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

modifyDesignMatrix <- function (rnaCond, riboCond)
{
    if (ncol(rnaCond) != ncol(riboCond))
        stop("rna-seq and ribo-seq must have the same number of conditions")
    numCond <- ncol(rnaCond)
    numRNASmps <- nrow(rnaCond)
    numRiboSmps <- nrow(riboCond)
    ### expand rna covariate vector with 0s
    rnaExpansion <- matrix(rep(rep(0,numCond), numRNASmps), nrow=numRNASmps)
    rnaCond <- cbind(rnaCond, as.data.frame(rnaExpansion))
    ### expand ribo covariate vector by repeating the same vector
    riboCond <- cbind(riboCond, riboCond)
    ### combine rna and ribo design matrix
    combinedCond <- rbind(rnaCond, riboCond)
    formula(combinedCond)
}

DESeq2Rex <- function (rnaCntTable, riboCntTable, rnaCond, riboCond)
{
    ## combine counts
    combCntTbl <- cbind(rnaCntTable, riboCntTable)

    if (ncol(rnaCond) != ncol(riboCond))
        stop("rna-seq and ribo-seq must have the same number of conditions")

    ### rnaCond <- cbind(intercept=factor(1), rnaCond)
    ### riboCond <- cbind(intercept=factor(1), riboCond)
    numCond <- ncol(rnaCond)
    numRNASmps <- nrow(rnaCond)
    numRiboSmps <- nrow(riboCond)
    ### expand rna covariate vector with 0s
    rnaExpansion <- matrix(factor(rep(rep(0,numCond), numRNASmps)), nrow=numRNASmps)
    rnaCond <- cbind(rnaCond, as.data.frame(rnaExpansion))
    numExtendedCols <- length(colnames(rnaCond))
    colnames(rnaCond)[(numCond+1):numExtendedCols] <- paste0("extra",seq(numCond))
    ### expand ribo covariate vector by repeating the same vector
    riboCond <- cbind(riboCond, riboCond)
    colnames(riboCond)[(numCond+1):numExtendedCols] <- paste0("extra",seq(numCond))
    ### combine rna and ribo design matrix
    combinedCond <- rbind(rnaCond, riboCond)
    extendedConds <- colnames(combinedCond)
    fmla <- as.formula(paste("~", paste(extendedConds, collapse= "+")))
    dds <- DESeqDataSetFromMatrix(countData = combCntTbl,
                                  colData = combinedCond,
                                  design = fmla)

    ## apply new design matrix with combine count table to DESeq2
    dds <- DESeq(dds)
    res <- results(dds)
    res
}

edgeRRex <- function (rnaCntTable, riboCntTable, rnaCond, riboCond)
{
    library(edgeR)

    ## combine design matrix
    ## order of data: control RNA samples, case RNA samples,
    ##                control Ribo samples, case Ribo samples
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps),
                          rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(rnaCntTable)), rep(1, ncol(riboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(rnaCntTable) + numCtlRiboSmps),
                         rep(1, numCaseRiboSmps)))

    ## combine counts
    combCntTbl <- cbind(rnaCntTable, riboCntTable)

    dge <- DGEList(counts = combCntTbl)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~condition+dataType+RiboDiff)
    dge <- estimateDisp(dge, design)

    ## glmFit and glmLRT
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)
    topGenes <- topTags(lrt, n=Inf)
}

voomRex <- function (rnaCntTable, riboCntTable, rnaCond, riboCond)
{
    library(edgeR)

    numCaseRNASmps <- ncol(rnaCntTable) - numCtlRNASmps

    group1 <- factor(c(rep(1,numCtlRNASmps), rep(2,numCaseRNASmps)))

    design1 <- model.matrix(~group1)

    numCaseRiboSmps <- ncol(riboCntTable) - numCtlRiboSmps

    group2 <- factor(c(rep(1,numCtlRiboSmps), rep(2,numCaseRiboSmps)))

    design2 <- model.matrix(~group2)

    ## combine design matrix
    ## order of data: control RNA samples, case RNA samples,
    ##                control Ribo samples, case Ribo samples
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps),
                          rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(RNACntTable)), rep(1, ncol(RiboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(RNACntTable) + numCtlRiboSmps),
                         rep(1, numCaseRiboSmps)))

    ## combine counts
    combCntTbl <- cbind(RNACntTable, RiboCntTable)

    dge <- DGEList(counts = combCntTbl)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~condition+dataType+RiboDiff)
    v <- voom(dge, design, plot=FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    topGenes <- topTable(fit, coef=ncol(design), number=Inf)
}

riborex <- function (rnaCntTable, riboCntTable, rnaCond, riboCond, engine)
{
    if(missing(engine)) {
      message("No engine selected, DESeq2 will be used by default")
    }
    if (engine == "DESeq2") {
        DESeq2Rex(rnaCntTable, riboCntTable, rnaCond, riboCond)
    }
    else if (engine == "edgeR") {
        edgeRRex(rnaCntTable, riboCntTable, rnaCond, riboCond)
    }
    else if (engine == "Voom") {
        voomRex(rnaCntTable, riboCntTable, rnaCond, riboCond)
    }
}
