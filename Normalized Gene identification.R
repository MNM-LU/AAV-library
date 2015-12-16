#suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(multicore))
library(plyr)
library(Hmisc)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
load.list <- read.table("loadlist.txt", header = FALSE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
colnames(load.list) <- c("Name", "BaseName","GroupName")
load.list <- rbind(load.list,c("completeLibraryRanges","","totalLib"))
load.list <- load.list[-grep("Untreat",load.list$Name),]

select.Cases <- c(unlist(sapply(load.list$Name, function(x) grep(x,in.names.all), simplify = TRUE)))

(in.names.all <- in.names.all[select.Cases])
grouping <- data.frame(Sample=gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.names.all)))),
                       Group=load.list[match(names(select.Cases),load.list$Name),"GroupName"],
                       stringsAsFactors = FALSE)

loadRDS <- function(in.name) {
  #in.name <- in.names.all[3]
  this.sample <- readRDS(in.name)
  this.name <- gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.name))))
  this.group <- grouping[match(this.name,grouping$Sample),"Group"]
  mcols(this.sample) <- cbind(mcols(this.sample),data.frame(Sample = this.name, Group=this.group,stringsAsFactors = FALSE))
  return(this.sample)
}

out.range <- lapply(in.names.all, loadRDS)
#do.call(sum,mcols(out.range)$tCount)
#out.range <- list(out.range[[1]][1:10000],out.range[[2]][1:10000])
readCounts <- lapply(out.range, function(x) sum(mcols(x)$RNAcount))
maxCount <- max(unlist(readCounts))
readCounts <- lapply(readCounts, function(x) maxCount/x)

makeNormCount <- function(inIndex){
  thisRange <- out.range[[inIndex]]
  mcols(thisRange)$RNAcount <- mcols(thisRange)$RNAcount*readCounts[[inIndex]]
  return(thisRange)
  }

out.range <- lapply(1:length(readCounts), makeNormCount)

out.range <- do.call(GAlignmentsList,unlist(out.range))
(out.range <- cbind(unlist(out.range))[[1]])
mcols(out.range)$NormCount <- 1

out.range.split <- split(out.range,c(mcols(out.range)$Group))
out.range.split <- lapply(out.range.split, function(x) split(x,seqnames(x)))
out.range.split <- mclapply(out.range.split, function(x) lapply(x, function(y) split(y,names(y))), mc.preschedule = TRUE, mc.cores = detectCores())

MergeCounts <- function(inRanges) {
outRanges <- inRanges[1]
mcols(outRanges) <- data.frame(structure=mcols(inRanges)$structure[1],
                               Group=mcols(inRanges)$Group[1],
                               LV=sum(mcols(inRanges)$LV*mcols(inRanges)$tCount)/sum(mcols(inRanges)$tCount),
                               mCount=sum(mcols(inRanges)$mCount),
                               tCount=sum(mcols(inRanges)$tCount),
                               BC=length(inRanges),
                               RNAcount=sum(mcols(inRanges)$RNAcount),
                               NormCount=log2(sum(mcols(inRanges)$RNAcount)+1)*length(inRanges),
                               stringsAsFactors=FALSE)
return(outRanges)
}

out.range.split <- lapply(out.range.split, function(x) mclapply(x, function(y) lapply(y,MergeCounts), mc.preschedule = TRUE, mc.cores = detectCores()))

out.range.split <- mclapply(out.range.split,function(x) unlist(do.call(GAlignmentsList,unlist(x)), use.names=FALSE), mc.preschedule = TRUE, mc.cores = detectCores())
out.range.split <- unlist(do.call(GAlignmentsList,unlist(out.range.split)), use.names=FALSE)
saveRDS(out.range.split, file="data/normalizedSampleRanges.RDS")
