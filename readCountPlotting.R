#suppressPackageStartupMessages(library(SummarizedExperiment))
load("alignedLibraries.rda")
load("LUTdna.rda")

makeTable <- function(in.range){
  table.analysis <- data.table(as.character(seqnames(in.range)), start(in.range)+(qwidth(in.range)/2), mcols(in.range)$RNAcount, 1L)
  setkey(table.analysis, V1) #Add V2 to allow for AA separation
  table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1)] #Add V2 to allow for AA separation
  #table.analysis.bin[,V2:=(V2+2)/3]
 return(table.analysis.bin)
}

load("completeLibraryRanges.rda")

table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L)
setkey(table.analysis, V1) #Add V2 to allow for AA separation        
table.analysis.bin <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1)] #Add V2 to allow for AA separation
#table.analysis.bin[,V2:=(V2+2)/3]
library.table <- table.analysis.bin

load("output/found.Nr08_100x_SN.rda")
Nr08_100x_SN.table <- makeTable(found.Nr08_100x_SN)


load("output/found.Nr08_100x_Ctx.rda")
Nr08_100x_Ctx.table <- makeTable(found.Nr08_100x_Ctx)


load("output/found.Nr08_100x_Str.rda")
Nr08_100x_Str.table <- makeTable(found.Nr08_100x_Str)

load("output/found.Nr08_100x_Thal.rda")
Nr08_100x_Thal.table <- makeTable(found.Nr08_100x_Thal)

load("output/found.Cells_293T_1000x.rda")
Cells_293T_1000x.table <- makeTable(found.Cells_293T_1000x)

save(library.table, Nr08_100x_SN.table, Nr08_100x_Ctx.table, Nr08_100x_Str.table, Nr08_100x_Thal.table, file="output/RNAtablesCompleteBin.rda")

complete.ranges.subset <- complete.ranges[1:100]

Union(granges(allFragments.ranges), complete.ranges.subset,
      ignore.strand=TRUE)
mapq_filter <- function(features, reads, ignore.strand, inter.feature) {
  require(GenomicAlignments) # needed for parallel evaluation
  Union(features, reads[mcols(reads)$mapq >= 20],
        ignore.strand=ignore.strand,
        inter.feature=inter.feature)
}

in.colData <- DataFrame("Sample1")
in.assay <- as.matrix(mcols(complete.ranges.subset)$tCount)

summedCount <- SummarizedExperiment(assays=SimpleList(counts=in.assay), rowRanges=granges(complete.ranges.subset),
                     colData=in.colData)
assays(summedCount)$counts



se <- summarizeOverlaps(features=granges(complete.ranges.subset), reads=complete.ranges.subset,
                        mode="IntersectionStrict",
                        singleEnd=TRUE,
                        ignore.strand=TRUE)

assays(se)$counts






se <- summarizeOverlaps(features=granges(allFragments.ranges), reads=complete.ranges,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE)
