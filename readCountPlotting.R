#suppressPackageStartupMessages(library(SummarizedExperiment))
load("alignedLibraries.rda")
load("LUTdna.rda")

makeTable <- function(in.range){
  table.analysis <- data.table(as.character(seqnames(in.range)), start(in.range)+(qwidth(in.range)/2), mcols(in.range)$RNAcount, 1L)
  setkey(table.analysis, V1, V2)        
  table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1,V2)]
  table.analysis.bin[,V2:=(V2+2)/3]
 return(table.analysis.bin)
}

load("completeLibraryRanges.rda")

table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L)
setkey(table.analysis, V1, V2)        
table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1,V2)]
table.analysis.bin[,V2:=(V2+2)/3]
library.table <- table.analysis.bin

load("output/found.Nr08_100x_SN.rda")
Nr08_100x_SN.table <- makeTable(found.Nr08_100x_SN)


load("output/found.Nr08_100x_Ctx.rda")
Nr08_100x_Ctx.table <- makeTable(found.Nr08_100x_Ctx)


load("output/found.Nr08_100x_Str.rda")
Nr08_100x_Str.table <- makeTable(found.Nr08_100x_Str)

load("output/found.Nr08_100x_Thal.rda")
Nr08_100x_Thal.table <- makeTable(found.Nr08_100x_Thal)

save(library.table, Nr08_100x_SN.table, Nr08_100x_Ctx.table, Nr08_100x_Str.table, Nr08_100x_Thal.table, file="output/RNAtables.rda")