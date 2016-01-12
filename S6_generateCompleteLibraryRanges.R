
#' ---
#' title: "Generate a complete library range object"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow clusters every read from each unique barcode from the library and matches them to the relevant ranges.  


suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(devtools))

load("data/alignedLibraries.rda")
load("data/LUTdna.rda")
load("data/multipleContfragmentsComplete.rda")
output.Table$fragment <- LUT.dna$Sequence[as.integer(output.Table$LUTnr)]
setkey(output.Table,fragment)

range.idx <- data.table(fragment=names(allFragments.ranges), 
                        idxFrag=1:length(allFragments.ranges), key="fragment")
output.Table <- output.Table[range.idx, nomatch=0, allow.cartesian=TRUE]

foundFragments.ranges <- allFragments.ranges[output.Table$idxFrag]
output.Table[,c("Reads","fragment","idxFrag"):=NULL]
mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges), output.Table)
  
saveRDS(foundFragments.ranges, file="output/completeLibraryRanges.rds")

devtools::session_info()