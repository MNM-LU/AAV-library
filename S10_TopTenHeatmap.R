
#' ---
#' title: "Top 10 heatmap analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This is the final script presenting top 10 candidates as heatmap plots.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
opts_chunk$set(fig.width = 5, fig.height = 9) #Full height 11
opts_chunk$set(comment = NA)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
library(grid)

#'Setup parameters
#'===================
all.samples <- readRDS("data/normalizedSampleRanges.RDS")
# fill.values <- c("CNS1000x_Ctx" = rgb(38,64,135, maxColorValue = 255), 
#                  "CNS1000x_Str" = rgb(157,190,217, maxColorValue = 255))
# filterBC <- FALSE
# filterAnimal <- FALSE
# BCadjustPlot <- FALSE
# NormalizePlot <- TRUE


#'Selection of relevant samples
#'===================
all.samples <- all.samples[-grep("4wks",mcols(all.samples)$Group)]
all.samples <- all.samples[-grep("PrimN_1000x_RNA",mcols(all.samples)$Group)]
select.samples <- all.samples[!(mcols(all.samples)$Group %in% "totalLib")]
trim.names <- data.table(GeneName=levels(seqnames(select.samples)))
trim.names[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
trim.names$GeneName <- gsub("/","-",trim.names$GeneName)
trim.names$GeneName <- gsub("_","-",trim.names$GeneName)
levels(seqnames(select.samples)) <- trim.names$GeneName
mcols(select.samples) <- cbind(mcols(select.samples), DataFrame(SeqCount=0))
select.samples.list <- split(select.samples, mcols(select.samples)$Group)
setCount <- function(inRange){
  #inRange <- select.samples.list[[3]]
  mcols(inRange[match(names(table(names(inRange))),names(inRange))])$SeqCount = as.integer(table(names(inRange)))
  inRange <- inRange[mcols(inRange)$SeqCount>=1]
    return(inRange)
}

select.samples.list <- lapply(select.samples.list,setCount)



select.samples <- unlist(do.call(GAlignmentsList,unlist(select.samples.list)), use.names=FALSE)


#'Selection of top ten fragments per sample
#'===================
o = order(-mcols(select.samples)$NormCount)
select.samples <- select.samples[o]
# sampleList <- unique(mcols(select.samples)$Group)
sampleList <- c("CNS100x_Str","CNS100x_Th","CNS100x_Ctx","CNS100x_SN")
#sampleList <- c("CNS1000x_Str","CNS1000x_Th","CNS1000x_Ctx","CNS1000x_SN")
#sampleList <- sampleList[c(4,3,1,2,5,7,6,9)]
#sampleList <- c("PrimN_100x","PrimN_1000x","293T_100x","293T_1000x")

topTen.list <- table(unlist(lapply(sampleList, function(x) names(select.samples[mcols(select.samples)$Group == x][1:min(10,length(select.samples[mcols(select.samples)$Group == x]))]))))
select.samples.top <- select.samples[names(select.samples) %in% names(topTen.list)]
full.coldata <- data.frame(Sample=sampleList, stringsAsFactors = TRUE)
row.names(full.coldata) <- sampleList
mcols(select.samples.top)$NormCount <- log2(mcols(select.samples.top)$RNAcount+1)




for (sample.name in sampleList) {
  #sample.name <- sampleList[1]
  this.sample <- select.samples.top[mcols(select.samples.top)$Group == sample.name]
  sub.matrix <- as.matrix(as.integer(mcols(this.sample)$NormCount, ncol = 1))
  row.names(sub.matrix) <- paste(seqnames(this.sample),(start(this.sample)+2)/3,mcols(this.sample)$structure, sep="_")
  colnames(sub.matrix) <- sample.name
  if (exists("full.matrix")) {
    full.matrix <- merge(full.matrix,sub.matrix,by="row.names",all.x=TRUE, all.y=TRUE, incomparables = NA)
    row.names(full.matrix) <- full.matrix$Row.names
    full.matrix$Row.names <- NULL
  } else {
    full.matrix <- sub.matrix
  }
}

full.matrix[is.na(full.matrix)] <- 0
#apply(full.matrix, 2, scale, center = FALSE) 
fullGene.analysis <- DESeqDataSetFromMatrix(countData = full.matrix,
                                            colData = full.coldata,
                                            design = ~ Sample)

norm.counts <- assay(fullGene.analysis)
pheatmap(norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)
pheatmap(norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)



# library("RColorBrewer")
# sampleDists <- dist(t(assay(fullGene.analysis.out)))
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- fullGene.analysis.out$Condition
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)


