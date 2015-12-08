#suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))

# in.BCs <- readFastq("data/barcodes_2015-11-05_AAVlibrary_complete.fastq.gz")
# BC.table <- rev(sort(table(sread(in.BCs))))
# BC.table <- BC.table/mean(BC.table)


load("completeLibraryRanges.rda")
libSum <- as.numeric(sum(mcols(complete.ranges)$tCount))
table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L, mcols(complete.ranges)$tCount)
setkey(table.analysis, V1) #Add V2 to allow for AA separation        
table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4), ReadNorm=sum(V5)), by=list(V1)] #Add V2 to allow for AA separation
table.analysis.bin$ReadNormZ <- scale(log2(table.analysis.bin$ReadNorm), center = TRUE)
#table.analysis.bin[,V2:=(V2+2)/3]
library.table <- table.analysis.bin



# load("alignedLibraries.rda")
# load("LUTdna.rda")


makeTable <- function(in.range){
  numRatio <- sum(mcols(in.range)$RNAcount)/libSum
  table.analysis <- data.table(as.character(seqnames(in.range)), start(in.range)+(qwidth(in.range)/2), mcols(in.range)$RNAcount, 1L, mcols(in.range)$RNAcount/numRatio)
  setkey(table.analysis, V1) #Add V2 to allow for AA separation        
  table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4), ReadNorm=sum(V5)), by=list(V1)] #Add V2 to allow for AA separation
  table.analysis.bin$ReadNormZ <- scale(log2(table.analysis.bin$ReadNorm), center = TRUE)
  #table.analysis.bin[,V2:=(V2+2)/3]
  return(table.analysis.bin)
}

in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
for (in.name in in.names.all){
assign(gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.name)))), makeTable(readRDS(in.name)))
}

out.names <- gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.names.all))))
out.names <- c(out.names,"library.table")
save(list = out.names, file="data/RNAtablesCompleteBin.rda")