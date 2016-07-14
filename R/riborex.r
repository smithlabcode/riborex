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

combineDesignMatrix <- function(rnaCond, riboCond) {

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (ncol(rnaCond) != ncol(riboCond))
      stop("rna- and ribo-seq data must have the same number of conditions")

  numCond <- ncol(rnaCond)
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)

  ### expand rna covariate vector with 0s
  rnaExpansion <- matrix(factor(rep(rep(0,numCond+1), numRNASmps)), nrow=numRNASmps)
  rnaCond <- cbind(rnaCond, as.data.frame(rnaExpansion))
  numExtendedCols <- length(colnames(rnaCond))
  colnames(rnaCond)[(numCond+1):numExtendedCols] <- paste0("extra",seq(numCond+1))

  ### expand ribo covariate vector by repeating the same vector
  riboCond <- cbind(riboCond, intercept=factor(1), riboCond)
  colnames(riboCond)[(numCond+1):numExtendedCols] <- paste0("extra",seq(numCond+1))

  ### combine rna and ribo design matrix
  combinedCond <- rbind(rnaCond, riboCond)

  extendedConds <- paste0("combinedCond$", colnames(combinedCond))
  fmla <- as.formula(paste("~", paste(extendedConds, collapse= "+")))
  model.matrix(fmla)
}

dataFrameToDesignMatrix <- function(cond) {
  if (!is.data.frame(cond)) cond <- data.frame(cond = cond)
  conditions <- paste0("cond$", colnames(cond))
  fmla <- as.formula(paste("~", paste(conditions, collapse= "+")))
  model.matrix(fmla)
}

DESeq2Rex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                      contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("rna- and ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (ncol(rnaCond) != ncol(riboCond))
    stop("rna- and ribo-seq data must have the same number of conditions")
  
  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  numCond <- ncol(rnaCond)
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)

  ### combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  ### expand rna covariate vector with 0s
  rnaExpansion <- matrix(factor(rep(rep(0,numCond+1), numRNASmps)),
                         nrow=numRNASmps)
  rnaCond <- cbind(rnaCond, as.data.frame(rnaExpansion))
  numExtendedCols <- length(colnames(rnaCond))
  colnames(rnaCond)[(numCond+1):numExtendedCols] <- paste0("extra",
                                                           seq(numCond+1))

  ### expand ribo covariate vector by repeating the same vector
  riboCond <- cbind(riboCond, intercept=factor(1), riboCond)
  colnames(riboCond)[(numCond+1):numExtendedCols] <- paste0("extra",
                                                            seq(numCond+1))

  ### combine rna and ribo design matrix
  combinedCond <- rbind(rnaCond, riboCond)
  extendedConds <- colnames(combinedCond)
  fmla <- as.formula(paste("~", paste(extendedConds, collapse= "+")))
  dds <- DESeqDataSetFromMatrix(countData = combCntTbl,
                                colData = combinedCond,
                                design = fmla)

  ## apply new design matrix with combined count table to DESeq2
  dds <- DESeq(dds)
  if(is.null(contrast)) {
    res <- results(dds)
  } else {
    res <- results(dds, contrast=contrast)
  }
  res
}

edgeRRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                     contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (ncol(rnaCond) != ncol(riboCond))
    stop("RNA- and Ribo-seq data must have the same number of conditions")
  
  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  ## combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  dge <- DGEList(counts = combCntTbl)
  dge <- calcNormFactors(dge)
  design <- combineDesignMatrix(rnaCond, riboCond)
  dge <- estimateDisp(dge, design)

  ## glmFit and glmLRT
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast=contrast)
  topGenes <- topTags(lrt, n=Inf)
  topGenes
}

edgeRDRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                      contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (ncol(rnaCond) != ncol(riboCond))
    stop("RNA- and Ribo-seq data must have the same number of conditions")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  ## estimate dispersion from RNA-seq data
  dge.rna <- DGEList(counts = rnaCntTable)
  dge.rna <- calcNormFactors(dge.rna)
  design.rna <- dataFrameToDesignMatrix(rnaCond)
  dge.rna <- estimateDisp(dge.rna, design.rna)

  ## estimate dispersion from Ribo-seq data
  dge.ribo <- DGEList(counts = riboCntTable)
  dge.ribo <- calcNormFactors(dge.ribo)
  design.ribo <- dataFrameToDesignMatrix(riboCond)
  dge.ribo <- estimateDisp(dge.ribo, design.ribo)

  ## combine dispersions
  dispersion.rna <- getDispersion(dge.rna)
  dispersion.ribo <- getDispersion(dge.ribo)
  dispersion <- matrix(c(rep(dispersion.rna, ncol(rnaCntTable)),
                         rep(dispersion.ribo, ncol(riboCntTable))),
                         nrow=dim(rnaCntTable)[1])
  ## combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)
  ## combine size factors
  combFactors <- c(dge.rna$samples$norm.factors,
                   dge.ribo$samples$norm.factors)
  ## create new DGE based on combined count table
  dge <- DGEList(counts = combCntTbl, norm.factors = combFactors)
  ## combine design matrix
  design <- combineDesignMatrix(rnaCond, riboCond)
  ## glmFit and glmLRT
  fit <- glmFit(dge, design=design, dispersion=dispersion)
  lrt <- glmLRT(fit, contrast=contrast)
  topGenes <- topTags(lrt, n=Inf)
  topGenes
}

voomRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                    contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (ncol(rnaCond) != ncol(riboCond))
    stop("RNA- and Ribo-seq data must have the same number of conditions")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) > minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) > minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  ## combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  dge <- DGEList(counts = combCntTbl)
  dge <- calcNormFactors(dge)
  design <- combineDesignMatrix(rnaCond, riboCond)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  if(!is.null(contrast)) {
    fit <- contrasts.fit(fit, contrasts = contrast)
  }
  fit <- eBayes(fit)

  topGenes <- topTable(fit, coef=ncol(design), number=Inf)
  topGenes
}

riborex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                    engine="DESeq2", contrast=NULL, minMeanCount=1) {

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
    stop ("Error: unrecognized engine name")
  }
}
