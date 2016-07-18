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

  message("combining design matrix")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (!identical(colnames(rnaCond), colnames(riboCond)))
      stop("RNA- and Ribo-seq data must have the same set of conditions")

  numCond <- ncol(rnaCond)
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)

  ### expand rna covariate vector with 0s
  expansion.rna <- rnaCond
  for(i in 1:ncol(expansion.rna)) {
    expansion.rna[,i] <- rnaCond[1,i]
  }
  expansion.rna <- cbind(0, as.data.frame(expansion.rna))
  rnaCond <- cbind(rnaCond, expansion.rna)
  colnames(rnaCond)[(numCond+1):ncol(rnaCond)] <- paste0("EXTRA",
                                                           seq(numCond+1))

  ### expand ribo covariate vector by repeating the same vector
  riboCond <- cbind(riboCond, 1, riboCond)
  colnames(riboCond)[(numCond+1):ncol(riboCond)] <- paste0("EXTRA",
                                                            seq(numCond+1))

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
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (!identical(colnames(rnaCond), colnames(riboCond)))
      stop("RNA- and Ribo-seq data must have the same set of conditions")

  if (ncol(rnaCntTable) != nrow(rnaCond))
    stop(paste("RNA-seq count table must have the",
               "same number of samples as in rnaCond"))

  if (ncol(riboCntTable) != nrow(riboCond))
    stop(paste("Ribo-seq count table must have the",
               "same number of samples as in riboCond"))

  if (minMeanCount < 1)
    stop("minMeanCount must at least be 1")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) >= minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) >= minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  numCond <- ncol(rnaCond)
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)

  ### combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  message("combining design matrix")

  ### expand rna covariate vector with 0s
  expansion.rna <- rnaCond
  for(i in 1:ncol(expansion.rna)) {
    expansion.rna[,i] <- rnaCond[1,i]
  }
  expansion.rna <- as.data.frame(cbind(0, expansion.rna))
  rnaCond <- cbind(rnaCond, expansion.rna)
  colnames(rnaCond)[(numCond+1):ncol(rnaCond)] <- paste0("EXTRA",
                                                           seq(numCond+1))

  ### expand ribo covariate vector by repeating the same vector
  riboCond <- cbind(riboCond, factor(1), riboCond)
  colnames(riboCond)[(numCond+1):ncol(riboCond)] <- paste0("EXTRA",
                                                            seq(numCond+1))

  ### combine rna and ribo design matrix
  combinedCond <- rbind(rnaCond, riboCond)
  extendedConds <- colnames(combinedCond)
  fmla <- as.formula(paste("~", paste(extendedConds, collapse= "+")))
  dds <- DESeqDataSetFromMatrix(countData = combCntTbl,
                                colData = combinedCond,
                                design = fmla)

  message("applying DESeq2 to modified design matrix")

  ## apply new design matrix with combined count table to DESeq2
  dds <- DESeq(dds)
  if(is.null(contrast)) {
    res <- results(dds)
  } else {
    res <- results(dds, contrast=contrast)
  }

  ## order results by gene names
  res <- res[order(rownames(res)),]
  res
}

edgeRRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                     contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (!identical(colnames(rnaCond), colnames(riboCond)))
      stop("RNA- and Ribo-seq data must have the same set of conditions")

  if (ncol(rnaCntTable) != nrow(rnaCond))
    stop(paste("RNA-seq count table must have the",
               "same number of samples as in rnaCond"))

  if (ncol(riboCntTable) != nrow(riboCond))
    stop(paste("Ribo-seq count table must have the",
               "same number of samples as in riboCond"))

  if (minMeanCount < 1)
    stop("minMeanCount must at least be 1")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) >= minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) >= minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  ## combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  dge <- DGEList(counts = combCntTbl)
  dge <- calcNormFactors(dge)
  design <- combineDesignMatrix(rnaCond, riboCond)
  dge <- estimateDisp(dge, design)

  message("applying edgeR to modified design matrix")

  ## glmFit and glmLRT
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast=contrast)
  topGenes <- topTags(lrt, n=Inf)

  ## order results by gene names
  topGenes <- topGenes[order(rownames(topGenes)),]
  topGenes
}

edgeRDRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                      contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (!identical(colnames(rnaCond), colnames(riboCond)))
      stop("RNA- and Ribo-seq data must have the same set of conditions")

  if (ncol(rnaCntTable) != nrow(rnaCond))
    stop(paste("RNA-seq count table must have the",
               "same number of samples as in rnaCond"))

  if (ncol(riboCntTable) != nrow(riboCond))
    stop(paste("Ribo-seq count table must have the",
               "same number of samples as in riboCond"))

  if (minMeanCount < 1)
    stop("minMeanCount must at least be 1")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) >= minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) >= minMeanCount)
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

  message("applying edgeR to modified design matrix")

  ## glmFit and glmLRT
  fit <- glmFit(dge, design=design, dispersion=dispersion)
  lrt <- glmLRT(fit, contrast=contrast)
  topGenes <- topTags(lrt, n=Inf)

  ## order results by gene names
  topGenes <- topGenes[order(rownames(topGenes)),]
  topGenes
}

voomRex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                    contrast=NULL, minMeanCount=1) {

  ### input validation
  if (!identical(rownames(rnaCntTable), rownames(riboCntTable)))
    stop ("RNA- and Ribo-seq data must have the same set of genes")

  if (!is.data.frame(rnaCond)) rnaCond <- data.frame(cond = rnaCond)
  if (!is.data.frame(riboCond)) riboCond <- data.frame(cond = riboCond)

  if (!identical(colnames(rnaCond), colnames(riboCond)))
      stop("RNA- and Ribo-seq data must have the same set of conditions")

  if (ncol(rnaCntTable) != nrow(rnaCond))
    stop(paste("RNA-seq count table must have the",
               "same number of samples as in rnaCond"))

  if (ncol(riboCntTable) != nrow(riboCond))
    stop(paste("Ribo-seq count table must have the",
               "same number of samples as in riboCond"))

  if (minMeanCount < 1)
    stop("minMeanCount must at least be 1")

  ### filter out low read count
  keep.rna <- which(rowMeans(rnaCntTable) >=  minMeanCount)
  keep.ribo <- which(rowMeans(riboCntTable) >=  minMeanCount)
  keep <- intersect(keep.rna, keep.ribo)
  rnaCntTable <- rnaCntTable[keep,]
  riboCntTable <- riboCntTable[keep,]

  ## combine counts
  combCntTbl <- cbind(rnaCntTable, riboCntTable)

  dge <- DGEList(counts = combCntTbl)
  dge <- calcNormFactors(dge)
  design <- combineDesignMatrix(rnaCond, riboCond)

  message("applying Voom to modified design matrix")

  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  if(!is.null(contrast)) {
    fit <- contrasts.fit(fit, contrasts = contrast)
  }
  fit <- eBayes(fit)

  topGenes <- topTable(fit, coef=ncol(design), number=Inf)

  ## order results by gene names
  topGenes <- topGenes[order(rownames(topGenes)),]
  topGenes
}

riborex <- function(rnaCntTable, riboCntTable, rnaCond, riboCond,
                    engine="DESeq2", contrast=NULL, minMeanCount=1) {

  if (engine == "DESeq2") {
    message("DESeq2 mode selected")
    DESeq2Rex(rnaCntTable, riboCntTable, rnaCond, riboCond,
              contrast, minMeanCount)
  }
  else if (engine == "edgeR") {
    message("edgeR mode selected")
    edgeRRex(rnaCntTable, riboCntTable, rnaCond, riboCond,
             contrast, minMeanCount)
  }
  else if (engine == "edgeRD") {
    message("edgeRD mode selected")
    edgeRDRex(rnaCntTable, riboCntTable, rnaCond, riboCond,
              contrast, minMeanCount)
  }
  else if (engine == "Voom") {
    message("Voom mode selected")
    voomRex(rnaCntTable, riboCntTable, rnaCond, riboCond,
            contrast, minMeanCount)
  }
  else {
    stop ("Error: unrecognized engine name")
  }
}
