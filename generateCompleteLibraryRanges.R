load("alignedLibraries.rda")
load("LUTdna.rda")
load("singleContfragments.rda")
output.Table.single <- output.Table
load("multipleContfragments.rda")
output.Table <- rbind(output.Table, output.Table.single)
rm(output.Table.single)
output.Table <- na.omit(output.Table)

output.Table$fragment <- LUT.dna$Sequence[as.integer(output.Table$LUTnr)]
#output.Table <- output.Table[1:1000,]
matchRange <- function(idxFrag) {
  #idxFrag <- 23
  machRanges <- which(names(allFragments.ranges) == output.Table$fragment[idxFrag])
  return(cbind(machRanges,idxFrag))
}
match.ranges.list <- mclapply(1:nrow(output.Table), matchRange, mc.preschedule = TRUE, mc.cores = detectCores())
match.ranges <- do.call(rbind, match.ranges.list)
complete.ranges <- allFragments.ranges[match.ranges[,1]]
mcols(complete.ranges) <- c(mcols(complete.ranges), output.Table[match.ranges[,2],2:6])
save(complete.ranges, file="completeLibraryRanges.rda")