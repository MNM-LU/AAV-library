suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
library("pheatmap")
attach("data/RNAtablesCompleteBin.rda")
sample.names <- ls("file:data/RNAtablesCompleteBin.rda")
sample.names <- c("library.table","total.infectiveLib",
      "RatNr1_100x_Str_7","RatNr7_100x_Str_3","RatNr8_100x_Str_11",
      "RatNr1_100x_Th_8","RatNr8_100x_Th_12","RatNr7_100x_Th_4",          
      "RatNr1_100x_Ctx_6","RatNr7_100x_Ctx_2","RatNr8_100x_Ctx_10",         
      "RatNr1_100x_SN_5","RatNr7_100x_SN_1","RatNr8_100x_SN_9",           
      "primNeuronsNr7_100x_cDNA_29","Cells293Nr3_100x_cDNA_27",
      "RatNr15_1000x_Str_15","RatNr19_1000x_Str_22","RatNr20_1000x_Str_24","RatNr21_1000x_Str_19",
      "RatNr15_1000x_Th_16","RatNr19_1000x_Th_23","RatNr20_1000x_Th_25","RatNr21_1000x_Th_20",
      "RatNr15_1000x_Ctx_14","RatNr19_1000x_Ctx_21","RatNr21_1000x_Ctx_18",
      "RatNr15_1000x_SN_13","RatNr21_1000x_SN_17",
      "primNeuronsNr6_1000x_cDNA_28","Cells293Nr2_1000x_cDNA_26")

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
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xStr",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xThal",
                                           "1000xCtx",
                                           "1000xCtx",
                                           "1000xCtx",
                                           "1000xSN",
                                           "1000xSN",
                                           "1000xPrimN",
                                           "1000x293T"), stringsAsFactors = TRUE)
row.names(full.coldata) <- sample.names
full.coldata$Tmp <- NA
full.coldata <- full.coldata[grepl("1000x",full.coldata$Condition),]
full.coldata$Tmp <- NULL
for (sample.name in row.names(full.coldata)) {
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
#rm(sub.matrix,this.sample,sample.name,sample.names)

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
#fullGene.analysis.out <- fullGene.analysis.out[,o]

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



