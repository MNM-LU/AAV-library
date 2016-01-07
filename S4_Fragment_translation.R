
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
suppressPackageStartupMessages(library(scales))

matchMethod="lv" #"hamming" "lv" "osa" "dl"

strt1<-Sys.time()

load("data/LUTdna.rda")

fragments.file <- "data/fragments_AAVlibrary_complete.fastq.gz"
barcodes.file <- "data/barcodes_AAVlibrary_complete.fastq.gz"

#' Make CustomArray reference index for Bowtie2
#' ============================
#+ Making Bowtie index.......


LUT.fa <- tempfile(pattern = "LUT_14aa_", tmpdir = tempdir(), fileext = ".fa")
LUT.seq = ShortRead(DNAStringSet(LUT.dna$Sequence), BStringSet(1:length(LUT.dna$Names)))
writeFasta(LUT.seq,LUT.fa)

LUT.14aa.fa <- tempfile(pattern = "LUT_14aa_", tmpdir = tempdir(), fileext = ".fa")
LUT.14aa.seq = ShortRead(DNAStringSet(LUT.dna[grep("14aa",LUT.dna$Structure),Sequence]), BStringSet(1:length(LUT.dna[grep("14aa",LUT.dna$Structure),Names])))
writeFasta(LUT.14aa.seq,LUT.14aa.fa)

LUT.22aa.fa <- tempfile(pattern = "LUT_22aa_", tmpdir = tempdir(), fileext = ".fa")
LUT.22aa.seq = ShortRead(DNAStringSet(LUT.dna[grep("22aa",LUT.dna$Structure),Sequence]), BStringSet(LUT.dna[grep("22aa",LUT.dna$Structure),Names]))
writeFasta(LUT.22aa.seq,LUT.22aa.fa)

# rm(LUT.14aa.seq,LUT.22aa.seq,LUT.seq)



#'Save unique fragments as fasta file
#'===================
reads.trim <- readFastq(fragments.file)
unique.reads <- unique(sread(reads.trim))

#unique.reads <- unique.reads[sample(length(unique.reads), 25000)]

unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))
fragments.unique.fa <- tempfile(pattern = "FragUnique_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(unique.reads,fragments.unique.fa)


#'Align against the library using blast
#'===================

blast.db <- tempfile(pattern = "blastDB_", tmpdir = tempdir(), fileext = ".db")
blast.out <- tempfile(pattern = "blastOut_", tmpdir = tempdir(), fileext = ".txt")

sys.out <-  system(paste("makeblastdb -in ", LUT.fa,
                         " -out ",blast.db," -dbtype nucl -title 'LUT' -parse_seqids 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("blastn database generation")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:(nrow(sys.out)),], format = "markdown")

sys.out <-  system(paste("export SHELL=/bin/sh; cat ",fragments.unique.fa," | parallel --block ",floor(length(unique.reads)/detectCores()),
                         " --recstart '>' --pipe blastn -max_target_seqs 10 -word_size 7",
                         " -num_threads 2 -outfmt 10 -db ", blast.db,
                         " -query - > ", blast.out, " 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) #-word_size 7

table.blastn <- data.table(read.table(blast.out, header = FALSE, skip = 0, sep=";",
                                     stringsAsFactors = FALSE, fill=FALSE),keep.rownames=FALSE)

if (length(grep("Warning",table.blastn$V1)) != 0) {
warnings.out <- unique(table.blastn[grep("Warning",table.blastn$V1),])
table.blastn <- table.blastn[-grep("Warning",table.blastn$V1),]
setnames(warnings.out,"V1", c("blastn Warnings"))
invisible(warnings.out[" "] <- " ")
knitr::kable(warnings.out[1:(nrow(warnings.out)),], format = "markdown")
}

table.blastn[,c("Reads","LUTnr","identity","alignmentLength","mismatches",
                 "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue","bitScore") := tstrsplit(V1,",",fixed=TRUE),]
table.blastn[,c("V1","identity","alignmentLength","gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue"):=NULL]

table.blastn[,Reads:= as.character(sread(unique.reads)[as.integer(Reads)])]
table.blastn[,bitScore:= as.numeric(bitScore)]
table.blastn[,mismatches:= as.numeric(mismatches)]
sort.order <- order(table.blastn$Reads,-table.blastn$bitScore)
table.blastn <- table.blastn[sort.order]
table.blastn.topHit <- table.blastn[!duplicated(table.blastn$Reads)]

reads.BC <- readFastq(barcodes.file)

full.table <- data.table(Reads=as.character(sread(reads.trim)),
                         BC=as.character(sread(reads.BC)),
                         key="Reads")

full.table <- merge(full.table,table.blastn.topHit, by="Reads", all.x = FALSE)

#' Starcode based barcode reduction
#' ============================
#+ Reducing barcodes.......

out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")

system(paste("gunzip -c ",barcodes.file," | starcode -t ",detectCores()-1," --print-clusters -d",
             1," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
       intern = TRUE, ignore.stdout = FALSE)

table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                                     stringsAsFactors = FALSE, fill=FALSE),keep.rownames=TRUE) #, nrows = 1000
table.BC.sc[,V2 := NULL]

table.BC.sc <- table.BC.sc[, strsplit(as.character(V3),",",fixed=TRUE), by=rn]

SC.droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
print(paste("Dropped BCs in Starcode:", SC.droppedBC))

#rm(reads.BC,reads.trim)

setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))

setkey(table.BC.sc,BC)

#' Replacing barcodes with Starcode reduced versions
#' ============================

setkey(full.table,BC)

full.table <- merge(full.table,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
rm(table.BC.sc)

setnames(full.table,c("BC","scBC"),c("oldBC","BC"))

setkey(full.table,BC)

RetainedBC <- length(unique(full.table$oldBC))
scBC <- length(unique(full.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))


table.frag <- data.table(as.data.frame((rev(sort(table(full.table$oldBC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
knitr::kable(table.frag, format = "markdown")

table.frag <- data.table(as.data.frame((rev(sort(table(full.table$BC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
knitr::kable(table.frag, format = "markdown")

invisible(full.table[,oldBC:=NULL])


#' Splitting reads into single-read and multi-read barcodes
#' ============================
#+ Splitting Reads.......

count.list <- table(full.table$BC)
full.table <- full.table[order(full.table$BC),]
temp.table.multi <- full.table[full.table$BC %in% names(count.list[count.list!=1]),]
temp.table.single <- full.table[full.table$BC %in% names(count.list[count.list==1]),]
temp.table.single[,c("mCount","tCount"):=1]
temp.table.single$Mode <- "Amb"
setkeyv(temp.table.multi,c("BC","LUTnr"))
key(temp.table.multi)

temp.table.multi[,mismatches:= as.numeric(mismatches)]
temp.table.multi[,c("bitScore","mismatches" ,"tCount"):= list(mean(bitScore),median(mismatches), .N), by=key(temp.table.multi)]
temp.table.multi$Mode <- "Def"
temp.table.multi <- unique(temp.table.multi)

print("Utilized reads.......")
print(nrow(full.table))
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
print(length(unique(temp.table.multi$BC)))

#' Calculate consensus alignment of chimeric barcodes
#' ============================
#+ Calculation consensus reads .....

table.blastn <- table.blastn[table.blastn$Reads %in% temp.table.multi$Reads,]

# calculate.consensus <- function(Reads,LUTnr,bitScore,tCount){
#   #invisible(lapply(c("Reads","LUTnr","bitScore","tCount"), function(x) assign(x, temp.table.multi[temp.table.multi$BC %in% "AAATCTGTGCGGGTTTAAGC",x, with=FALSE], envir = .GlobalEnv)))
#     
#   group.table <- data.table(Reads,tCount,key="Reads")
#   group.table[, c("mCount","tCount"):=list(tCount,sum(tCount))]
#   score.table <- table.blastn[table.blastn$Reads %in% group.table$Reads,]
#   group.table <- merge(group.table,score.table,by="Reads")
#   group.table[,c("bitScore","mismatches" ,"mCount"):= list(max(bitScore),median(mismatches), sum(mCount)), by="LUTnr"]
#   group.table <- group.table[group.table$mCount == max(group.table$mCount),]
#   group.table <- group.table[which.max(group.table$bitScore),]
#   
#   if (max(group.table$mCount) == 1){
#     group.table$Mode <- "Amb"
#   } else {
#     group.table$Mode <- "Def"
#   }
# 
#   return(group.table[1,])
# }
# #temp.table.multi <- temp.table.multi[1:10000,]
# #temp.table.multi.sub <- temp.table.multi[4100000:4100100,]
# strt8 <- Sys.time()
# temp.table.multi.sub[,c("Reads","tCount","mCount","LUTnr","mismatches","bitScore","Mode"):=calculate.consensus(Reads,LUTnr,bitScore,tCount), by="BC"]
# print(Sys.time()-strt8)


#Second alternative

# calculate.consensus <- function(group.table){
#   #group.table <- temp.table.multi[temp.table.multi$BC %in% "GTTTGCGTGTTGCCGCGTGC",]
#   group.table[, c("mCount","tCount"):=list(tCount,sum(tCount))]
#   score.table <- table.blastn[table.blastn$Reads %in% group.table$Reads,]
#   group.table <- merge(group.table,score.table,by="Reads", allow.cartesian=TRUE)
#   group.table[,c("bitScore","mismatches" ,"mCount"):= list(max(bitScore),median(mismatches), sum(mCount)), by="LUTnr"]
#   group.table <- group.table[group.table$mCount == max(group.table$mCount),]
#   group.table <- group.table[which.max(group.table$bitScore),]
#   
#   if (max(group.table$mCount) == 1){
#     group.table$Mode <- "Amb"
#   } else {
#     group.table$Mode <- "Def"
#   }
#   
#   return(group.table[1,])
# }
# #temp.table.multi <- temp.table.multi[4100000:4100100,]
# setkey(temp.table.multi,BC)
# 
# temp.table.multi[,c("LUTnr","bitScore","mismatches"):=NULL]
# temp.table.multi.table <- split(temp.table.multi,temp.table.multi$BC)
# temp.table.multi.table <- mclapply(temp.table.multi.table, calculate.consensus, 
#                                        mc.preschedule = TRUE, mc.cores = detectCores())
# temp.table.multi <- rbindlist(temp.table.multi.table)

#Third alternative



calculate.consensus <- function(Reads,LUTnr,bitScore,tCount){
  #invisible(lapply(c("Reads","LUTnr","bitScore","tCount"), function(x) assign(x, temp.table.multi[temp.table.multi$BC %in% "AAATCTGTGCGGGTTTAAGC",x, with=FALSE], envir = .GlobalEnv)))
    
  group.table <- data.table(Reads,tCount,key="Reads")
  group.table[, c("mCount","tCount"):=list(tCount,sum(tCount))]
  score.table <- table.blastn[table.blastn$Reads %in% group.table$Reads,]
  group.table <- merge(group.table,score.table,by="Reads")
  group.table[,c("bitScore","mismatches" ,"mCount"):= list(max(bitScore),median(mismatches), sum(mCount)), by="LUTnr"]
  group.table <- group.table[group.table$mCount == max(group.table$mCount),]
  group.table <- group.table[which.max(group.table$bitScore),]
  
  if (max(group.table$mCount) == 1){
    group.table$Mode <- "Amb"
  } else {
    group.table$Mode <- "Def"
  }

  return(group.table[1,])
}

#temp.table.multi <- temp.table.multi[1:10000,]
#temp.table.multi.sub <- temp.table.multi[4100000:4100100,]
split.table <- data.table(BC=unique(temp.table.multi$BC), key="BC")
split.table$Cut <- cut(1:nrow(split.table),detectCores()/2)
temp.table.multi <- merge(temp.table.multi,split.table, by="BC")
setkey(temp.table.multi,"Cut")

rm(table.blastn.topHit, split.table,reads.BC,reads.trim, unique.reads,full.table,LUT.dna)
gc()

print(system.time(temp.table.multi.table <- split(temp.table.multi,temp.table.multi$Cut)))

strt8 <- Sys.time()
# temp.table.multi.table <- mclapply(temp.table.multi.table, function(x) x[,c("Reads","tCount","mCount","LUTnr","mismatches","bitScore","Mode"):=calculate.consensus(Reads,LUTnr,bitScore,tCount), by="BC"], 
#                                    mc.preschedule = FALSE, mc.cores = detectCores()/2)

temp.table.multi.table <- mclapply(temp.table.multi.table, function(x) x[,j=list(calculate.consensus(Reads,LUTnr,bitScore,tCount)), by="BC"], 
                                   mc.preschedule = FALSE, mc.cores = detectCores()/2)

temp.table.multi <- rbindlist(temp.table.multi.table)
temp.table.multi[,"Cut":=NULL]
print(Sys.time()-strt8)



setkeyv(temp.table.multi,c("BC","LUTnr"))
#temp.table.multi <- unique(temp.table.multi)
temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean)

print(paste("Total number of definitive Barcodes:", length(grep("Def",temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", length(grep("Amb",temp.table.multi.consensus$Mode))))
print(paste("Total number of single-read Barcodes:", nrow(temp.table.single)))

output.Table <- rbind(temp.table.multi.consensus,temp.table.single)
save(output.Table, file="data/multipleContfragmentsNew.rda")

print("Total analysis time:")
print(Sys.time()-strt1)
devtools::session_info()