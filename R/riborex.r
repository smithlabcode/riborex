DESeq2Rex <- function (rnaCntTable, riboCntTable, rnaCondition, riboCondition)
{
    library(DESeq2)

    # combine counts
    combCntTbl <- cbind(rnaCntTable, riboCntTable)

    # combine design matrix
    # order of data: control RNA samples, case RNA samples,
    #                control Ribo samples, case Ribo samples
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps),
                   rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(RNACntTable)), rep(1, ncol(RiboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(RNACntTable) + numCtlRiboSmps),
                  rep(1, numCaseRiboSmps)))

    # combine design matrix
    combColData <- data.frame(condition = condition)
    combColData$dataType <- dataType
    combColData$RiboDiff <- RiboDiff
    rownames(combColData) <- colnames(combCntTbl)

    dds <- DESeqDataSetFromMatrix(countData = combCntTbl,
                                  colData = combColData,
                                  design = ~condition + dataType + RiboDiff)

    # apply new design matrix with combine count table to DESeq2
    dds <- DESeq(dds)
    res <- results(dds)
    res
}

riborex <- function (rnaCntTable, riboCntTable, rnaCondition, riboCondition, method)
{
    if (method == "DESeq2") {
        DESeq2Rex(rnaCntTable, riboCntTable, rnaCondition, riboCondition)
    }
}
