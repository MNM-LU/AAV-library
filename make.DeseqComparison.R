suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
library("pheatmap")
attach("data/RNAtablesCompleteBin.rda")
sample.names <- ls("file:data/RNAtablesCompleteBin.rda")
o = c(3,23,9,17,21,10,18,22,7,15,19,8,16,20,6,2,13,14,11,12,5,1)
sample.names <- sample.names[o]
full.coldata <- data.frame(Condition=rbind("PlasmidLib",
                                           "InfectiveAAV",
                                           "100xStr",
                                           "100xStr",
                                           "100xStr",
                                           "100xThal",
                                           "100xThal",
                                           "100xThal",
                                           "100xCtx",
                                           "100xCtx",
                                           "100xCtx",
                                           "100xSN",
                                           "100xSN",
                                           "100xSN",
                                           "100xPrimN",
                                           "100x293T",
                                           "1000xStr",
                                           "1000xThal",
                                           "1000xCtx",
                                           "1000xSN",
                                           "1000xPrimN",
                                           "1000x293T"), stringsAsFactors = TRUE)
row.names(full.coldata) <- sample.names
for (sample.name in sample.names) {
  #sample.name=sample.names[3]
  this.sample <- get(sample.name)
  sub.matrix <- as.matrix(this.sample$ReadCount, ncol = 1)
  row.names(sub.matrix) <- this.sample$V1
  colnames(sub.matrix) <- sample.name
  if (exists("full.matrix")) {
    full.matrix <- merge(full.matrix,sub.matrix,by="row.names",all.x=TRUE, all.y=TRUE, incomparables = NA)
    row.names(full.matrix) <- full.matrix$Row.names
    full.matrix$Row.names <- NULL
  } else {
    full.matrix <- sub.matrix
  }
}
detach("file:data/RNAtablesCompleteBin.rda")
rm(sub.matrix,this.sample,sample.name,sample.names,o)

full.matrix[is.na(full.matrix)] <- 0
gene.names <- data.table(rownames(full.matrix))
gene.names[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(V1, ",", fixed=TRUE)]
rownames(full.matrix) <- gene.names$GeneName

fullGene.analysis <- DESeqDataSetFromMatrix(countData = full.matrix,
                              colData = full.coldata,
                              design = ~ Condition)

fullGene.analysis <- collapseReplicates(fullGene.analysis, full.coldata$Condition, renameCols = TRUE)

fullGene.analysis.log2 <- normTransform(fullGene.analysis) # defaults to log2(x+1)
# fullGene.analysis.vst <- varianceStabilizingTransformation(fullGene.analysis, blind=FALSE)
# fullGene.analysis.rlog <- rlog(fullGene.analysis, blind=FALSE)

log2.norm.counts <- assay(fullGene.analysis.log2)
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE)
