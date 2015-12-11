suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
library("pheatmap")
attach("data/RNAtablesCompleteBin.rda")
sample.names <- ls("file:data/RNAtablesCompleteBin.rda")
o = c(3,53,11,47,51,34,13,48,52,35,7,45,49,32,9,46,50,33,43,44,6,2,23,30,36,40,17,25,31,37,41,18,19,29,38,15,21,27,28,4,1)
sample.names <- sample.names[o]
full.coldata <- data.frame(Condition=rbind("PlasmidLib",
                                           "InfectiveAAV",
                                           "100xStr",
                                           "100xStr",
                                           "100xStr",
                                           "100xStr",
                                           "100xThal",
                                           "100xThal",
                                           "100xThal",
                                           "100xThal",
                                           "100xCtx",
                                           "100xCtx",
                                           "100xCtx",
                                           "100xCtx",
                                           "100xSN",
                                           "100xSN",
                                           "100xSN",
                                           "100xSN",
                                           "100xMu",
                                           "100xSc",
                                           "100xPrimN",
                                           "100x293T",
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xCtx",
                                           "1000xCtx",
                                           "1000xCtx",
                                           "1000xCtx",
                                           "1000xSN",
                                           "1000xMu",
                                           "1000xSc",
                                           "1000xPrimN",
                                           "1000x293T"), stringsAsFactors = TRUE)
row.names(full.coldata) <- sample.names
for (sample.name in sample.names) {
  #sample.name=sample.names[3]
  this.sample <- get(sample.name)
  sub.matrix <- as.matrix(as.integer(this.sample$ReadNorm, ncol = 1))
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
full.matrix <- full.matrix+1
gene.names <- data.table(rownames(full.matrix))
gene.names[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(V1, ",", fixed=TRUE)]
rownames(full.matrix) <- gene.names$GeneName

fullGene.analysis <- DESeqDataSetFromMatrix(countData = full.matrix,
                              colData = full.coldata,
                              design = ~ Condition)

fullGene.analysis <- collapseReplicates(fullGene.analysis, full.coldata$Condition, renameCols = TRUE)
fullGene.analysis.log2 <- normTransform(fullGene.analysis, f = log2, pc = 0)
#fullGene.analysis.vst <- varianceStabilizingTransformation(fullGene.analysis, blind=FALSE)
#fullGene.analysis.rlog <- rlog(fullGene.analysis, blind=FALSE)
fullGene.analysis.out <- fullGene.analysis.log2
 # defaults to log2(x+1)
# fullGene.analysis.vst <- varianceStabilizingTransformation(fullGene.analysis, blind=FALSE)


#o = c(14,13,11,12,8,10,9,7,5,6,2,4,3,1)
#o = c(14,13,11,12,8,10,9,7)
#o = c(5,6,2,4,3,1)
# o = c(11,12,8,10,9,7)
# fullGene.analysis.out <- fullGene.analysis.out[,o]

norm.counts <- assay(fullGene.analysis.out)
pheatmap(norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE)


library("RColorBrewer")
sampleDists <- dist(t(assay(fullGene.analysis.out)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- fullGene.analysis.out$Condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



