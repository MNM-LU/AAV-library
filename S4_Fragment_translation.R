
#' ---
#' title: "Barcoded plasmid library translation"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow clusters every read from each unique barcode and determines the consensus fragmnt from the CustomArray and the barcode fidelity. i.e., if the barcode was monoclonal or not..  


suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringdist))

matchMethod="lv" #"hamming" "lv" "osa" "dl"

strt1<-Sys.time()

load("data/LUTdna.rda")
reads.trim <- readFastq("data/fragments_2015-11-05_AAVlibrary_complete.fastq.gz")
reads.BC <- readFastq("data/barcodes_2015-11-05_AAVlibrary_complete.fastq.gz")

reads.table <- sread(reads.trim)
names(reads.table) <- sread(reads.BC)



#' Starcode based barcode reduction
#' ============================
#+ Reducing barcodes.......

in.name.BC.star <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")

writeFastq(reads.BC,in.name.BC.star,compress=TRUE)

system(paste("gunzip -c ",in.name.BC.star," | starcode -t ",detectCores()/2," --print-clusters -d",
             1," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
       intern = TRUE, ignore.stdout = FALSE)

table.BC.sc <- read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                          stringsAsFactors = FALSE, fill=FALSE) #, nrows = 1000
table.BC.sc$V2 <- NULL

#list.BC.sc <- split(table.BC.sc, rownames(table.BC.sc))
list.BC.sc.list <- mclapply(1:nrow(table.BC.sc), function(x) strsplit(as.character(table.BC.sc$V3[x]), ","), mc.preschedule = TRUE, mc.cores = detectCores(), mc.cleanup = TRUE)
names(list.BC.sc.list) <- rownames(table.BC.sc)
list.BC.sc.list <- lapply(list.BC.sc.list, rbind)
table.BC.sc <- data.table(cbind(unlist(list.BC.sc.list, use.names = TRUE)), keep.rownames=TRUE)
invisible(table.BC.sc[,rn:=gsub("[0-9]","",rn)])
SC.droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
print(paste("Dropped BCs in Starcode:", SC.droppedBC))

rm(reads.BC,reads.trim)

setnames(table.BC.sc,'V1','BC')
setnames(table.BC.sc,'rn','scBC')

setkey(table.BC.sc,BC)
table.BC.sc<- unique(table.BC.sc)
BCs.sc.counts <- table(table.BC.sc$BC)
BCs.sc.counts.single <- BCs.sc.counts[BCs.sc.counts == 1]
table.BC.sc <- table.BC.sc[table.BC.sc$BC %in% names(BCs.sc.counts.single)]

temp.table <- data.table(BC=names(reads.table),Reads=as.character(reads.table))
setkey(temp.table,BC)
temp.table <- merge(temp.table,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
rm(table.BC.sc)

setnames(temp.table,'BC','oldBC')
setnames(temp.table,'scBC','BC')
setkey(temp.table,BC)

RetainedBC <- length(unique(temp.table$oldBC))
scBC <- length(unique(temp.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))


table.frag <- data.table(as.data.frame((rev(sort(table(temp.table$oldBC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
knitr::kable(table.frag, format = "markdown")

table.frag <- data.table(as.data.frame((rev(sort(table(temp.table$BC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
knitr::kable(table.frag, format = "markdown")

invisible(temp.table[,oldBC:=NULL])
reads.table <- DNAStringSet(temp.table$Reads)
names(reads.table) <- temp.table$BC
#rm(temp.table)



strt2<-Sys.time()

#' Aligning unique fragments to the CustomArray reference
#' ============================
#+ Aligning to reference.......


LUT.fa <- tempfile(pattern = "LUT_", tmpdir = tempdir(), fileext = ".fa")
LUT.seq = ShortRead(DNAStringSet(LUT.dna$Sequence), BStringSet(1:length(LUT.dna)))
writeFasta(LUT.seq,LUT.fa)

bowtieIDX <- tempfile(pattern = "IDX_LUT_", tmpdir = tempdir(), fileext = "")



sample(nrow(temp.table.small), 1000)
temp.table.small <- data.table(Reads=unique(temp.table$Reads), key="Reads")
#temp.table.small <- temp.table.small[11000:12000,]
# temp.table.small.tmp <- temp.table.small
# temp.table.small.tmp$LUTnr <- NULL
# 

temp.table.small <- temp.table.small[order(temp.table.small$Reads)]
temp.table.small$LUTnr <- match(temp.table.small$Reads,LUT.dna$Sequence)
temp.table.small.exact <- temp.table.small[!is.na(temp.table.small$LUTnr),]
strt5<-Sys.time()
temp.table.small$LUTnr <- amatch(temp.table.small$Reads, LUT.dna$Sequence, method=matchMethod, 
       maxDist = 30, matchNA = FALSE, useBytes = TRUE, nthread = getOption("sd_num_thread"))
print(Sys.time()-strt5)

strt3<-Sys.time()
temp.table.small$LUTseq <- LUT.dna$Sequence[temp.table.small$LUTnr]
temp.table.small <- temp.table.small[!is.na(temp.table.small$LUTseq)]
temp.table.small[,LV:= stringdist(Reads,LUTseq, method=matchMethod, nthread = getOption("sd_num_thread"))]

setkey(temp.table,Reads)
temp.table.out <- merge(temp.table,temp.table.small, by="Reads", all = FALSE, all.x = FALSE)
setkeyv(temp.table.out,c("BC","LUTnr"))
temp.table.out[,c("Reads","LUTseq"):=NULL]

#' Splitting reads into single-read and multi-read barcodes
#' ============================
#+ Splitting Reads.......

count.list <- table(temp.table.out$BC)
temp.table.multi <- temp.table.out[temp.table.out$BC %in% names(count.list[count.list!=1])]
temp.table.single <- temp.table.out[temp.table.out$BC %in% names(count.list[count.list==1])]
temp.table.single[,c("mCount","tCount"):=1]
key(temp.table.multi)

temp.table.multi[,c("LV","tCount"):= list(mean(LV), .N), by=key(temp.table.multi)]
temp.table.multi <- unique(temp.table.multi)

print("Utilized Barcodes.......")
print(nrow(temp.table.out))
print("Whereof single reads.......")
print(nrow(temp.table.single))

#' Splitting multi-read barcodes into clean and chimeric
#' ============================
#+ Splitting Clean Reads.......

setkeyv(temp.table.multi,"BC")
count.list <- table(temp.table.multi$BC)

temp.table.multi.clean <- temp.table.multi[temp.table.multi$BC %in% names(count.list[count.list==1])]
temp.table.multi <- temp.table.multi[temp.table.multi$BC %in% names(count.list[count.list!=1])]
temp.table.multi.clean[,mCount:=tCount]

print("Clean multi-read barcodes.......")
print(nrow(temp.table.multi.clean))
print("Chimeric multi-read barcodes.......")
print(nrow(temp.table.multi))

#' Calculate consensus alignment of chimeric barcodes
#' ============================
#+ Calculation consensus reads .....

calculate.consensus <- function(LUTnr,LV,tCount){
  group.table <- data.table(LUTnr,LV,tCount,key="LUTnr")
  group.table[, mCount:=tCount]
  group.table[, tCount:=sum(tCount)]
  if (max(tCount) == 1){
    group.table <- group.table[which.min(group.table$LV),]
  } else {
    group.table <- group.table[which.max(group.table$mCount),]
  }
  if (nrow(group.table) > 1){
    group.table[,LUTnr:=NA]
  }
  return(group.table[1,])
}


temp.table.multi[,c("LUTnr","LV","tCount","mCount"):=calculate.consensus(LUTnr,LV,tCount), by="BC"]
setkeyv(temp.table.multi,c("BC","LUTnr"))
temp.table.multi <- unique(temp.table.multi)
temp.table.multi <- temp.table.multi[!is.na(temp.table.multi$LUTnr)]
temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean)

output.Table <- temp.table.multi.consensus
save(output.Table, file="data/multipleContfragmentsNew.rda")
output.Table <- temp.table.multi.clean
save(output.Table, file="data/singleContfragmentsNew.rda")

print("Post-matching analysis time:")
print(Sys.time()-strt3)


print("Fragment translation time:")
print(Sys.time()-strt2)


print("Total analysis time:")
print(Sys.time()-strt1)
devtools::session_info()