
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

fragments.file <- "data/fragments_2015-11-05_AAVlibrary_subset.fastq.gz"
barcodes.file <- "data/barcodes_2015-11-05_AAVlibrary_subset.fastq.gz"

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
unique.reads = ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))
fragments.unique.fa <- tempfile(pattern = "FragUnique_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(unique.reads,fragments.unique.fa)


#'Align against the 14aa library using usearch
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
warnings.out <- unique(table.blastn[grep("Warning",table.blastn$V1),])
setnames(warnings.out,"V1", c("blastn Warnings"))
invisible(warnings.out[" "] <- " ")
knitr::kable(warnings.out[1:(nrow(warnings.out)),], format = "markdown")


table.blastn <- table.blastn[-grep("Warning",table.blastn$V1),]

table.blastn[,c("query id","subject id","% identity","alignment length","mismatches",
                 "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue","bit score") := tstrsplit(V1,",",fixed=TRUE),]
table.blastn[,V1:=NULL]

#'Align against the 14aa library using bowtie2
#'===================

bowtieIDX <- tempfile(pattern = "IDX_LUT_", tmpdir = tempdir(), fileext = "")
sys.out <-  system(paste("bowtie2-build",LUT.14aa.fa,bowtieIDX,  sep = " "), 
                   intern = TRUE, ignore.stdout = FALSE) 
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 build index of 14aa CustomArray fragments")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:30,], format = "markdown")

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --threads ",detectCores(),
                         " --very-sensitive",
                         " -x ", bowtieIDX, " -U ",fragments.file," -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to CustomArray fragments")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:(nrow(sys.out)),], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

fragment.ranges.14aa <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)

#'Align against the 22aa library using bowtie2
#'===================
bowtieIDX <- tempfile(pattern = "IDX_LUT_", tmpdir = tempdir(), fileext = "")
sys.out <-  system(paste("bowtie2-build",LUT.22aa.fa,bowtieIDX,  sep = " "), 
                   intern = TRUE, ignore.stdout = FALSE) 
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 build index of 22aa CustomArray fragments")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:30,], format = "markdown")

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --threads ",detectCores(),
                         " --very-sensitive",
                         " -x ", bowtieIDX, " -U ",fragments.file," -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to CustomArray fragments")
invisible(sys.out[" "] <- " ")
knitr::kable(sys.out[1:(nrow(sys.out)),], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

fragment.ranges.22aa <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)


#rm(LUT.14aa.fa,LUT.22aa.fa,bowtieIDX,name.bowtie,sys.out)
#'Merge reads and alignments
#'===================


reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)

full.table <- data.table(Reads=as.character(sread(reads.trim)),
                         BC=as.character(sread(reads.BC)),
                         ID=unlist(lapply(strsplit(as.character(ShortRead::id(reads.trim)), " "),"[",1)),
                         key="ID")

found.order <- as.integer(match(seqnames(fragment.ranges.14aa),LUT.dna$Names))
aligned.table.14aa <- data.table(ID=as.character(names(fragment.ranges.14aa)),
                                 LUTnr14=found.order,
                                 LUTseq14=as.character(LUT.dna$Sequence[found.order]),
                                 Cigar14=as.character(cigar(fragment.ranges.14aa)),
                                 key="ID")

found.order <- as.integer(match(seqnames(fragment.ranges.22aa),LUT.dna$Names))
aligned.table.22aa <- data.table(ID=as.character(names(fragment.ranges.22aa)),
                                 LUTnr22=found.order,
                                 LUTseq22=as.character(LUT.dna$Sequence[found.order]),
                                 Cigar22=as.character(cigar(fragment.ranges.22aa)),
                                 key="ID")

full.table.matched <- merge(full.table,aligned.table.14aa, by="ID", all.x = TRUE)
full.table.matched <- merge(full.table.matched,aligned.table.22aa, by="ID", all.x = TRUE)
full.table.matched[,ID := NULL]

found.order <- is.na(full.table.matched$LUTnr14) & is.na(full.table.matched$LUTnr22)
full.table.nonFound <- full.table.matched[found.order,]
full.table.matched <- full.table.matched[!found.order,]

print(paste("Percent of all reads successfully aligned:", percent(nrow(full.table.matched)/nrow(full.table))))

#rm(aligned.table.14aa,aligned.table.22aa,full.table,fragment.ranges.14aa,fragment.ranges.22aa,found.order)

reads.nonFound <- data.table(Reads=unique(full.table.nonFound$Reads), key="Reads")
#reads.nonFound <- reads.nonFound[sample(nrow(reads.nonFound), 1000),]
strt3<-Sys.time()
reads.nonFound$LUTnr <- amatch(reads.nonFound$Reads, LUT.dna$Sequence, method=matchMethod, 
                               maxDist = 8, matchNA = FALSE, useBytes = TRUE, nthread = getOption("sd_num_thread"))
print(Sys.time()-strt3)
reads.nonFound[,LUTseq := LUT.dna$Sequence[LUTnr]]
full.table.nonFound <- merge(full.table.nonFound,reads.nonFound, by="Reads", all=FALSE, all.x = FALSE)

#' Aligning unique fragments to the CustomArray reference
#' ============================
#+ Aligning to reference.......


temp.table.small <- full.table.matched #[sample(nrow(full.table.matched), 10000),]

strt3<-Sys.time()
temp.table.small[,LV14:= stringdist(Reads,LUTseq14, method=matchMethod, nthread = getOption("sd_num_thread"))]
temp.table.small[,LV22:= stringdist(Reads,LUTseq22, method=matchMethod, nthread = getOption("sd_num_thread"))]
full.table.nonFound[,LV:= stringdist(Reads,LUTseq, method=matchMethod, nthread = getOption("sd_num_thread"))]
temp.table.small.22aa <- temp.table.small[(temp.table.small$LV22 <= temp.table.small$LV14) | is.na(temp.table.small$LV14),]
temp.table.small.14aa <- temp.table.small[(temp.table.small$LV14 < temp.table.small$LV22) | is.na(temp.table.small$LV22),]
temp.table.small.22aa[,c("LV14","LUTseq14","LUTnr14","Cigar14") := NULL]
temp.table.small.14aa[,c("LV22","LUTseq22","LUTnr22","Cigar22") := NULL]
full.table.nonFound[,c("LV14","LUTseq14","LUTnr14","Cigar14","LV22","LUTseq22","LUTnr22") := NULL]
temp.table.small.22aa$LUTlib <- "22aa"
temp.table.small.14aa$LUTlib <- "14aa"
full.table.nonFound$LUTlib <- "Manual"
setnames(temp.table.small.22aa,c("LV22","LUTseq22","LUTnr22","Cigar22"),c("LV","LUTseq","LUTnr","Cigar"))
setnames(temp.table.small.14aa,c("LV14","LUTseq14","LUTnr14","Cigar14"),c("LV","LUTseq","LUTnr","Cigar"))
setnames(full.table.nonFound,c("Cigar22"),c("Cigar"))

temp.table.small <- rbind(temp.table.small.22aa,temp.table.small.14aa,full.table.nonFound)


setkeyv(temp.table.small,c("BC","LUTnr"))




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

rm(reads.BC,reads.trim)

setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))

setkey(table.BC.sc,BC)

#' Replacing barcodes with Starcode reduced versions
#' ============================

setkey(temp.table.small,BC)

temp.table <- merge(temp.table.small,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
rm(table.BC.sc)

setnames(temp.table,c("BC","scBC"),c("oldBC","BC"))

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


#' Splitting reads into single-read and multi-read barcodes
#' ============================
#+ Splitting Reads.......
temp.table.out <- temp.table
count.list <- table(temp.table.out$BC)
temp.table.multi <- temp.table.out[temp.table.out$BC %in% names(count.list[count.list!=1])]
temp.table.multi <- temp.table.multi[order(BC)]
temp.table.single <- temp.table.out[temp.table.out$BC %in% names(count.list[count.list==1])]
temp.table.single[,c("mCount","tCount"):=1]
setkeyv(temp.table.multi,c("BC","LUTnr"))
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
#temp.table.multi <- temp.table.multi[1:10000,]

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