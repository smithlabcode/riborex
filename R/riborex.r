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

edgeRRex <- function (rnaCntTable, riboCntTable, rnaCondition, riboCondition)
{
    library(edgeR)
    
    # combine design matrix
    # order of data: control RNA samples, case RNA samples,
    #                control Ribo samples, case Ribo samples
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps),
    rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(RNACntTable)), rep(1, ncol(RiboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(RNACntTable) + numCtlRiboSmps),
    rep(1, numCaseRiboSmps)))
    
    # combine counts
    combCntTbl <- cbind(RNACntTable, RiboCntTable)
    
    dge <- DGEList(counts = combCntTbl)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~condition+dataType+RiboDiff)
    dge <- estimateDisp(dge, design)
    
    # glmFit and glmLRT
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)
    topGenes <- topTags(lrt, n=Inf)
}

voomRex <- function (rnaCntTable, riboCntTable, rnaCondition, riboCondition)
{
    library(edgeR)
    # read in RNA Count Table
    RNACntTable <- read.table(RNACntFile, row.names = 1, header = TRUE)
    
    numCaseRNASmps <- ncol(RNACntTable) - numCtlRNASmps
    
    group1 <- factor(c(rep(1,numCtlRNASmps), rep(2,numCaseRNASmps)))
    
    design1 <- model.matrix(~group1)
    
    # read in Ribo Count Table
    RiboCntTable <- read.table(RiboCntFile, row.names = 1, header = TRUE)
    
    numCaseRiboSmps <- ncol(RiboCntTable) - numCtlRiboSmps
    
    group2 <- factor(c(rep(1,numCtlRiboSmps), rep(2,numCaseRiboSmps)))
    
    design2 <- model.matrix(~group2)
    
    # combine design matrix
    # order of data: control RNA samples, case RNA samples,
    #                control Ribo samples, case Ribo samples
    condition <- factor(c(rep(0, numCtlRNASmps), rep(1, numCaseRNASmps),
    rep(0, numCtlRiboSmps), rep(1, numCaseRiboSmps)))
    dataType <- factor(c(rep(0, ncol(RNACntTable)), rep(1, ncol(RiboCntTable))))
    RiboDiff <- factor(c(rep(0, ncol(RNACntTable) + numCtlRiboSmps),
    rep(1, numCaseRiboSmps)))
    
    # combine counts
    combCntTbl <- cbind(RNACntTable, RiboCntTable)
    
    dge <- DGEList(counts = combCntTbl)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~condition+dataType+RiboDiff)
    v <- voom(dge, design, plot=FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    topGenes <- topTable(fit, coef=ncol(design), number=Inf)
}

riborex <- function (rnaCntTable, riboCntTable, rnaCondition, riboCondition, method)
{
    if (method == "DESeq2") {
        DESeq2Rex(rnaCntTable, riboCntTable, rnaCondition, riboCondition)
    }
    else if (method == "edgeR") {
        edgeRRex(rnaCntTable, riboCntTable, rnaCondition, riboCondition)
    }
    else if (method == "voom") {
        voomRex(rnaCntTable, riboCntTable, rnaCondition, riboCondition)
    }
}
