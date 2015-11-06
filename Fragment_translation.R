suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(beanplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales)) #Gives the log2 ability to ggplot2
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(matrixStats))


LUT.dna <- read.table("Complete fragment list for Custom array 2015-02-10.txt", header = TRUE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
LUT.dna <- data.table(LUT.dna)
invisible(LUT.dna[,Sequence:=gsub("aacctccagagaggcaac","",Sequence)])
invisible(LUT.dna[,Sequence:=gsub("cagacaagcagctaccgca","",Sequence)])
invisible(LUT.dna[,Sequence:=toupper(Sequence)])


reads.trim <- readFastq("data/fragments_2015-11-05_AAVlibrary_complete.fastq.gz")
reads.BC <- readFastq("data/barcodes_2015-11-05_AAVlibrary_complete.fastq.gz")


# reads.BC.file <- "data/barcodes_2015-11-04_AAVlibrary.fastq.gz"
# 
# 
# setkey(barcodeTable,BC)
# 
# #' Starcode based barcode reduction
# #' ============================
# #+ Reducing barcodes.......
# 
# out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")
# system(paste("gunzip -c ",reads.BC.file," | starcode -t ",detectCores()/2," --print-clusters -d",
#              0," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
#        intern = TRUE, ignore.stdout = FALSE)
# 
# table.BC.sc <- read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
#                           stringsAsFactors = FALSE, fill=FALSE)
# table.BC.sc$V2 <- NULL
# list.BC.sc <- split(table.BC.sc, rownames(table.BC.sc))
# list.BC.sc <- lapply(list.BC.sc, function(x) strsplit(as.character(x), ","))
# list.BC.sc <- lapply(list.BC.sc, rbind)
# table.BC.sc <- data.table(cbind(unlist(list.BC.sc, use.names = TRUE)), keep.rownames=TRUE)
# invisible(table.BC.sc[,rn:=gsub("[0-9]","",rn)])
# output.table$droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
# print(paste("Dropped BCs in Starcode:", output.table$droppedBC[1]))
# 
# setnames(table.BC.sc,'V1','BC')
# setnames(table.BC.sc,'rn','scBC')
# 
# setkey(table.BC.sc,BC)
# table.BC.sc<- unique(table.BC.sc)
# BCs.sc.counts <- table(table.BC.sc$BC)
# BCs.sc.counts.single <- BCs.sc.counts[BCs.sc.counts == 1]
# table.BC.sc <- table.BC.sc[table.BC.sc$BC %in% names(BCs.sc.counts.single)]
# 
# barcodeTable <- merge(barcodeTable,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
# rm(table.BC.sc)
# rm(reads.BC)
# setnames(barcodeTable,'BC','oldBC')
# setnames(barcodeTable,'scBC','BC')
# setkey(barcodeTable,BC)
# 
# output.table$RetainedBC <- length(unique(barcodeTable$oldBC))
# output.table$scBC <- length(unique(barcodeTable$BC))
# print(paste("Original unique barcodes:", output.table$RetainedBC[1]))
# print(paste("SC reduced unique barcodes:", output.table$scBC[1]))
# 
# 
# table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$oldBC))))[1:10]), keep.rownames=TRUE)
# setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
# knitr::kable(table.frag, format = "markdown")
# 
# table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$BC))))[1:10]), keep.rownames=TRUE)
# setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
# knitr::kable(table.frag, format = "markdown")
# 
# invisible(barcodeTable[,oldBC:=NULL])

#reads.BC <- reads.BC[width(sread(reads.trim)) < 75L & width(sread(reads.trim)) > 38L]
reads.trim <- reads.trim[width(sread(reads.trim)) < 78L & width(sread(reads.trim)) > 38L]

source("retrieveFASTAQID.R")
FastQ1ID <- retrieveFASTAQID(reads.BC, PE=TRUE)
FastQ2ID <- retrieveFASTAQID(reads.trim, PE=TRUE)

hits <- intersect(FastQ2ID,FastQ1ID)

reads.BC <- reads.BC[match(hits,FastQ1ID)]
reads.trim <- reads.trim[match(hits,FastQ2ID)]
reads.table <- sread(reads.trim)
names(reads.table) <- sread(reads.BC)

BCs.counts <- table(names(reads.table))
BCs.single <- BCs.counts[BCs.counts == 1]
print("Utilized Barcodes.......")
print(length(BCs.counts))
print("Whereof single reads.......")
print(length(BCs.single))
reads.table.single <- reads.table[names(reads.table) %in%  names(BCs.single)]
reads.table.multiple <- reads.table[!(names(reads.table) %in%  names(BCs.single))]

reads.table.list <- split(reads.table.multiple, names(reads.table.multiple))
rm(reads.BC,reads.trim,FastQ1ID,FastQ2ID,hits,reads.table)

#reads.table.list <- reads.table.list[1:100]

#reads.table.list <- reads.table.list[mcmapply(length,reads.table.list,mc.preschedule = TRUE, mc.cores = detectCores()/2L) > 3]

#reads.table.list <- reads.table.list[sapply(reads.table.list,length) > 3]

#reads <- reads.table.list[[99]]

match.pair <- function(reads){
uniques <- rev(sort(table(reads)))
read.match <- adist(names(uniques),LUT.dna$Sequence)
read.match <- do.call(rbind,lapply(1:nrow(read.match), function(i) c(which.min(read.match[i,]), min(read.match[i,]))))
uniques <- data.table(cbind(read.match, uniques))
uniques <- split(uniques, uniques$V1)
uniques <- data.table(do.call(rbind, lapply(uniques,function(x) c((sum(x$V2*x$uniques)/sum(x$uniques)), sum(x$uniques)))), keep.rownames=TRUE)
uniques <- uniques[which.max(uniques$V2),]
uniques$V3 <- length(reads)
return(uniques)
}


#reads.sub <- reads.table.list[90:100]

output.Table = data.frame(matrix(vector(), length(reads.table.list), 5,
                                 dimnames=list(c(), c("LUTnr","LV","mCount","tCount", "BC"))),
                          stringsAsFactors=F)


for (startPos in seq(1,length(reads.table.list),1000)) {
reads.sub <- reads.table.list[startPos:min((startPos+999),length(reads.table.list))]
strt4<-Sys.time()
match.out <- mclapply(reads.sub, match.pair, mc.preschedule = TRUE, mc.cores = detectCores()-1L )
table.sub <- data.table(do.call(rbind,match.out))
table.sub$BC <- names(reads.sub)
firstRow <- sum(!is.na(output.Table[,1]))+1
output.Table[firstRow:(firstRow+nrow(table.sub)-1),] <- table.sub
save(output.Table, file="multipleContfragments.rda")
print("Total fragment translation time:")
print(Sys.time()-strt4)
}



match.pair.single <- function(read){
  read.match <- adist(read,LUT.dna$Sequence)
  min.LV <- min(read.match)
  read.match <- which.min(read.match)
  return(data.table(cbind(read.match, min.LV, 1L, 1L)))
}

reads.single.sub <- reads.table.single[300:310]

strt4<-Sys.time()
match.out.single <- mclapply(reads.single.sub, match.pair.single, mc.preschedule = TRUE, mc.cores = detectCores()/2L )
table.sub.single <- data.table(do.call(rbind,match.out.single))
table.sub.single$BC <- names(reads.single.sub)

print("Total fragment translation time:")
print(Sys.time()-strt4)

colnames(table.sub) <- c("LUTnr","LV","mCount","tCount", "BC")
colnames(table.sub.single) <- c("LUTnr","LV","mCount","tCount", "BC")

table.complete <- rbind(table.sub,table.sub.single)
setkey(table.complete, "BC")

table.complete