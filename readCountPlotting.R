#suppressPackageStartupMessages(library(SummarizedExperiment))
load("alignedLibraries.rda")
load("LUTdna.rda")

makeTable <- function(in.range){
  table.analysis <- data.table(as.character(seqnames(in.range)), start(in.range)+(qwidth(in.range)/2), mcols(in.range)$RNAcount, 1L)
  setkey(table.analysis, V1, V2) #Add V2 to allow for AA separation
  table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1, V2)] #Add V2 to allow for AA separation
  table.analysis.bin$normCount <- table.analysis.bin$ReadCount/mean(table.analysis.bin$ReadCount)
  table.analysis.bin[,V2:=(V2+2)/3]
 return(table.analysis.bin)
}

load("completeLibraryRanges.rda")

table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L)
setkey(table.analysis, V1, V2) #Add V2 to allow for AA separation        
table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1, V2)] #Add V2 to allow for AA separation
table.analysis.bin$normCount <- table.analysis.bin$ReadCount/mean(table.analysis.bin$ReadCount)
table.analysis.bin[,V2:=(V2+2)/3]
library.table <- table.analysis.bin

in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
for (in.name in in.names.all){
assign(gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.name)))), makeTable(readRDS(in.name)))
}

out.names <- gsub("-","_",gsub("found.","",gsub("(output/)", "", gsub("(.rds)", "", in.names.all))))
out.names <- c(out.names,"library.table")
save(list = out.names, file="data/RNAtables.rda")