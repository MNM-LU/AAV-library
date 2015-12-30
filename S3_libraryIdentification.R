
#' ---
#' title: "Library identification"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow aligns the library fragments to the full reference sequences using Bowtie2. 
suppressPackageStartupMessages(library(knitr)) 
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))


opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(comment = NA)

#'Load sequences
#'===================
LUT.dna <- read.table("data/SortedFragments_all.txt",
                      header = TRUE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
LUT.dna <- data.table(LUT.dna)

#'Remove constitutive backbone sequences
#'===================
invisible(LUT.dna[,Sequence:=gsub("aacctccagagaggcaacg","",Sequence)])
invisible(LUT.dna[,Sequence:=gsub("cagacaagcagctaccgca","",Sequence)])
invisible(LUT.dna[,Sequence:=toupper(Sequence)])
setkey(LUT.dna, "Sequence")
LUT.dna <- unique(LUT.dna)
LUT.dna$Names <- LUT.dna$Sequence
save(LUT.dna,file = "data/LUTdna.rda")

#'Split sequences based on linker and length 
#'===================

#output.Table$LUTseq <- LUT.dna$Sequence[as.numeric(output.Table$LUTnr)]
LUT.14aaG4S <- LUT.dna[substr(LUT.dna$Sequence,1,14) == "GAGGCGGAGGAAGT"]
LUT.remaining <- LUT.dna[!(substr(LUT.dna$Sequence,1,14) == "GAGGCGGAGGAAGT")]
LUT.14aaA5 <- LUT.remaining[substr(LUT.remaining$Sequence,1,14) == "CTGCTGCAGCAGCC"]
LUT.remaining <- LUT.remaining[!(substr(LUT.remaining$Sequence,1,14) == "CTGCTGCAGCAGCC")]
LUT.22aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 70L & 
                            substr(LUT.remaining$Sequence,1,2) == "CT"]
LUT.remaining <- LUT.remaining[!(nchar(LUT.remaining$Sequence) == 70L & 
                                   substr(LUT.remaining$Sequence,1,2) == "CT")]
LUT.14aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 46L & 
                            substr(LUT.remaining$Sequence,1,2) == "CT"]

#'Trim sequences
#'===================
LUT.14aa$Sequence <- substr(LUT.14aa$Sequence,3,44)
LUT.14aaG4S$Sequence <- substr(LUT.14aaG4S$Sequence,15,56)
LUT.14aaA5$Sequence <- substr(LUT.14aaA5$Sequence,15,56)
LUT.22aa$Sequence <- substr(LUT.22aa$Sequence,3,68)

#'Save fasta files for Bowtie alignments
#'===================

LUT.14aa.fa <- tempfile(pattern = "LUT_14aa_", tmpdir = tempdir(), fileext = "fa")
LUT.14aa.seq = ShortRead(DNAStringSet(LUT.14aa$Sequence), BStringSet(LUT.14aa$Names))
writeFasta(LUT.14aa.seq,LUT.14aa.fa)

LUT.14aaG4S.fa <- tempfile(pattern = "LUT_14aaG4s_", tmpdir = tempdir(), fileext = "fa")
LUT.14aaG4S.seq = ShortRead(DNAStringSet(LUT.14aaG4S$Sequence), BStringSet(LUT.14aaG4S$Names))
writeFasta(LUT.14aaG4S.seq,LUT.14aaG4S.fa)

LUT.14aaA5.fa <- tempfile(pattern = "LUT_14aaA5_", tmpdir = tempdir(), fileext = "fa")
LUT.14aaA5.seq = ShortRead(DNAStringSet(LUT.14aaA5$Sequence), BStringSet(LUT.14aaA5$Names))
writeFasta(LUT.14aaA5.seq,LUT.14aaA5.fa)

LUT.22aa.fa <- tempfile(pattern = "LUT_14aaA5_", tmpdir = tempdir(), fileext = "fa")
LUT.22aa.seq = ShortRead(DNAStringSet(LUT.22aa$Sequence), BStringSet(LUT.22aa$Names))
writeFasta(LUT.22aa.seq,LUT.22aa.fa)


#'Build Bowtie index
#'===================

seqs.original <- readFasta("input/DNA Libraries for Retrograde Transport_fPolChikv.fasta")

seqs.AA <- Biostrings::translate(sread(seqs.original), genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

source("functions/AAtoDNA.R")
seqs.optimized = ShortRead(DNAStringSet(sapply(seqs.AA, function(x) AAtoDNA(x, species="hsa"))), BStringSet(gsub("([ ])", "_", ShortRead::id(seqs.original))))

bowtie.fasta <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = ".fa")

writeFasta(seqs.optimized,bowtie.fasta)

bowtie.idx <- tempfile(pattern = "IDX_bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2-build",bowtie.fasta,bowtie.idx, "2>&1",  sep = " "), 
                   intern = TRUE, ignore.stdout = FALSE) 



#' Align fragments to reference
#' ============================
#+ Align to reference...

#' Align 14aa sequences
name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores(),
                         " --local --score-min 'C,0,-1' -f -a",
                         " -x ", bowtie.idx, " -U ", LUT.14aa.fa, " -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), 
                   intern = TRUE, ignore.stdout = FALSE) 


sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aa.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aa.ranges))
length(unique(names(frag14aa.ranges)))
length(unique(LUT.14aa$Sequence))

#' Align 14aaG4S sequences

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores(),
                         " --local --score-min 'C,0,-1' -f -a",
                         " -x ", bowtie.idx, " -U ", LUT.14aaG4S.fa, " -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE)

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aaG4S.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aaG4S.ranges))
length(unique(names(frag14aaG4S.ranges)))
length(unique(LUT.14aaG4S$Sequence))

#' Align 14aaA5 sequences

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores(),
                         " --local --score-min 'C,0,-1' -f -a",
                         " -x ", bowtie.idx, " -U ", LUT.14aaA5.fa, " -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aaA5.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aaA5.ranges))
length(unique(names(frag14aaA5.ranges)))
length(unique(LUT.14aaA5$Sequence))

#' Align 22aa sequences

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores(),
                         " --local --score-min 'C,0,-1' -f -a",
                         " -x ", bowtie.idx, " -U ",LUT.22aa.fa," -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -Sb ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag22aa.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag22aa.ranges))
length(unique(names(frag22aa.ranges)))
length(unique(LUT.22aa$Sequence))

#' Merge and annotate aligned sequences
#' ============================

mcols(frag14aa.ranges)$structure <- "14aa"
mcols(frag22aa.ranges)$structure <- "22aa"
mcols(frag14aaA5.ranges)$structure <- "14aaA5"
mcols(frag14aaG4S.ranges)$structure <- "14aaG4S"
allFragments.ranges <- append(frag14aa.ranges,frag22aa.ranges)
allFragments.ranges <- append(allFragments.ranges,frag14aaA5.ranges)
allFragments.ranges <- append(allFragments.ranges,frag14aaG4S.ranges)

save(allFragments.ranges, file="data/alignedLibraries.rda")

devtools::session_info()