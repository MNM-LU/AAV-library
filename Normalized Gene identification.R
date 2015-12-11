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
select.Cases <- c(grep("100x_Str",in.names.all),grep("100x_Th",in.names.all),grep("100x_Ctx",in.names.all),grep("100x_SN",in.names.all),
                  grep("100x_Mu",in.names.all),grep("100x_Sc",in.names.all),
                  grep("primNeuronsNr7_100x",in.names.all),grep("Cells293Nr3_100x",in.names.all),
                  grep("1000x_Str",in.names.all),grep("1000x_Th",in.names.all),grep("1000x_Ctx",in.names.all),grep("1000x_SN",in.names.all),
                  grep("100x_Mu",in.names.all),grep("100x_Sc",in.names.all),
                  grep("primNeuronsNr6_1000x_cDNA",in.names.all),grep("Cells293Nr2_1000x",in.names.all))
(in.names.all <- in.names.all[select.Cases])
in.names.all <- in.names.all[-c(1,6,11,16,26,32,37,43)]
grouping <- data.frame(Files=gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.names.all)))),
                       Group=rbind(rep.row("CNS100x_Str",4),
                                   rep.row("CNS100x_Th",4),
                                   rep.row("CNS100x_Ctx",4),
                                   rep.row("CNS100x_SN",4),
                                   "PerN100x_Mu",
                                   "PerN100x_SC",
                                   "PrimN_100x",
                                   "293T_100x",
                                   rep.row("CNS1000x_Str",5),
                                   rep.row("CNS1000x_Th",5),
                                   rep.row("CNS1000x_Ctx",4),
                                   rep.row("CNS1000x_SN",3),
                                   "PerN1000x_Mu",
                                   "PerN1000x_SC",
                                   "PrimN_1000x",
                                   "293T_1000x"))

loadRDS <- function(in.name) {
  this.sample <- readRDS(in.name)
  mcols(this.sample) <- cbind(mcols(this.sample),data.frame(Sample = gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.name)))),stringsAsFactors = FALSE))
  return(this.sample)
}

out.range <- lapply(in.names.all, loadRDS)
do.call(sum,mcols(out.range)$tCount)


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
mcols(out.range)$Sample <- sedit(mcols(out.range)$Sample,as.character(grouping$Files), as.character(grouping$Group), wild.literal=FALSE)
mcols(out.range)$NormCount <- 1

out.range.split <- split(out.range,c(mcols(out.range)$Sample))
out.range.split <- lapply(out.range.split, function(x) split(x,seqnames(x)))
out.range.split <- mclapply(out.range.split, function(x) lapply(x, function(y) split(y,names(y))), mc.cores = detectCores())

MergeCounts <- function(inRanges) {
outRanges <- inRanges[1]
mcols(outRanges)$NormCount <- log2(sum(mcols(inRanges)$RNAcount)+1)*length(inRanges)
return(outRanges)
}

out.range.split <- mclapply(out.range.split, function(x) lapply(x, function(y) lapply(y,MergeCounts)), mc.cores = detectCores())

out.range.split <- mclapply(out.range.split,function(x) unlist(do.call(GAlignmentsList,unlist(x)), use.names=FALSE), mc.cores = detectCores())
out.range.split <- unlist(do.call(GAlignmentsList,unlist(out.range.split)), use.names=FALSE)




makeTable <- function(in.range){
  numRatio <- sum(mcols(in.range)$RNAcount)/libSum
  table.analysis <- data.table(as.character(seqnames(in.range)), start(in.range)+(qwidth(in.range)/2), mcols(in.range)$RNAcount, 1L, mcols(in.range)$RNAcount/numRatio)
  setkey(table.analysis, V1) #Add V2 to allow for AA separation        
  table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4), ReadNorm=sum(V5)), by=list(V1)] #Add V2 to allow for AA separation
  table.analysis.bin$ReadNormZ <- scale(log2(table.analysis.bin$ReadNorm), center = TRUE)
  #table.analysis.bin[,V2:=(V2+2)/3]
  return(table.analysis.bin)
}






assign(gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.name)))), makeTable(readRDS(in.name)))


out.names <- gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.names.all))))
out.names <- c(out.names,"library.table")
save(list = out.names, file="data/RNAtablesCompleteBin.rda")