
#' ---
#' title: "Library analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This workflow brings together FastQ files containing barcodes and 5'/3' ends of a suitable insert and alignmen them using Bowtie2. It also includes starcode based false barcode reduction and a MapReduce based hierarchical clustering  
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



#+ setup, include=FALSE
opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(comment = NA)


config <- read.table("config.txt", header = FALSE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
colnames(config) <- c("Parameter", "Value")
#setwd("~/Dropbox (Bjorklund Lab)/Shared/NGS data")
script.dir <- normalizePath("../SharedFunctions/")
# source(file.path(script.dir, "buildFragLengthTable.R"))
# source(file.path(script.dir, "makePEfastq.R"))
# source(file.path(script.dir, "retrieveFASTAQID.R"))
# source(file.path(script.dir, "buildFragTableSC3.R"))

output.table <- data.frame(meanCount=numeric(),
                           sdCount=numeric(), 
                           meanP5=numeric(),
                           sdP5=numeric(),
                           meanP7=numeric(),
                           sdP7=numeric(),
                           meanWidth=numeric(),
                           sdWidth=numeric(),
                           q0=numeric(),
                           q25=numeric(),
                           q50=numeric(),
                           q75=numeric(),
                           q100=numeric(),
                           Name=character(),
                           SC=character(),
                           Reads=numeric(),
                           OrigBC=numeric(),
                           RetainedBC=numeric(),
                           scBC=numeric(),
                           droppedBC=numeric(),
                           uniqueFragFwd=numeric(),
                           uniqueFragRev=numeric(),
                           stringsAsFactors=FALSE) 

output.table[1,1:as.integer(ncol(output.table))] <- NA

#'Sequencing files
#'===================
knitr::kable(config, format = "markdown")
dataDir <- config$Value[1]
in.name.P5 <- file.path(dataDir, config$Value[2])
in.name.P7 <- file.path(dataDir, config$Value[3])
name.out <- config$Value[4]
paired.alignment <- as.logical(config$Value[5])

#'Analysis parameters
#'===================
bb.dir <- config$Value[6]
fragmentTemplate  <- config$Value[7]
output.table$SC <- config$Value[8]
run.subset <- as.logical(config$Value[9])
align.p7 <- as.logical(config$Value[10])
max.cores <- as.integer(config$Value[11])
subset.count <- as.integer(config$Value[12])

#'Script execution
#'===================
strt<-Sys.time()

id.backbone.L <- file.path(bb.dir, "Ltrim.fa") 
id.backbone.R <- file.path(bb.dir, "Rtrim.fa")
id.BC.L <- file.path(bb.dir, "BC-L.fa")
id.BC.R <- file.path(bb.dir, "BC-R.fa")
id.uncut <- file.path(bb.dir, "uncut.fa")


#' Selection of real amplicons
#' ============================
#+ Selecting real amplicons.......

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("-Xmx12g overwrite=true k=15 rcomp=f skipr2=t qhdist=0 maskmiddle=f hammingdistance=2 findbestmatch=f ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "GTATGTTGTTCTGGAGCGGGAGGGTGCTATTTTGCCTAGCGATAA", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t , ACAAGCAGCTACCGCAGATGTCAACACA           
# postLoxP on P5: GTATGTTGTTCTGGAGCGGGAGGGTGCTATTTTGCCTAGCGATAAGCTGATGTAGCC
# GFP from P7: CCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
# Cap from P7: AGACAAGCAGCTACCGCAGATGTCAACACACAAGGCGTTCTTCCAGGCATGGTCTGG

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)


colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7



#' Extraction of a subset
#' ============================
#+ Extracting subset.......

if (run.subset){
  suppressWarnings(sampler <- FastqSampler(gsub("([\\])", "", in.name.P5), subset.count, readerBlockSize=1e9, ordered = TRUE)) 
  set.seed(123); tmp.P5 <- yield(sampler)
  in.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
  writeFastq(tmp.P5,in.name.P5, compress=TRUE)
  rm(tmp.P5)
  suppressWarnings(sampler <- FastqSampler(gsub("([\\])", "", in.name.P7), subset.count, readerBlockSize=1e9, ordered = TRUE)) 
  set.seed(123); tmp.P7 <- yield(sampler)
  in.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
  writeFastq(tmp.P7,in.name.P7, compress=TRUE)
  rm(tmp.P7)
}

output.table$Reads <- as.integer(system(paste("gunzip -c ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE)) #Stores the read count utilized
print(paste("Utilized sequences:", output.table$Reads[1]))


#' Extraction of barcodes
#' ============================
#+ Extracting barcodes.......

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("-Xmx12g overwrite=true k=10 rcomp=f skipr2=t qhdist=0 maskmiddle=t hammingdistance=1 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "ATAACTTCGTATAATGTATGC", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t ,

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)


colnames(sys.out) <- c("bbduk2 Identification of real barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7


out.name.P5 <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=15 mink=15 hammingdistance=1 findbestmatch=t ",
                        "rcomp=f findbestmatch=f qhdist=0 minavgquality=0 maxns=0 minlength=18 ",
                        "maxlength=22 threads=", detectCores()," in=", shQuote(in.name.P5), 
                        " out=", out.name.P5," lliteral=", "GGCCTAGCGGCCGCTTTACTT",
                        " rliteral=", "ATAACTTCGTATAATGTATGC",
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) #" fliteral=",id.uncut,
sys.out <- as.data.frame(sys.out)

in.name.P5 <- out.name.P5


colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")
rm(sys.out)

reads.BC <- readFastq(in.name.P5)
sread(reads.BC)
output.table$OrigBC <- length(unique(sread(reads.BC)))
unique(sread(reads.BC))
barcodeTable <- data.table(ID=as.character(ShortRead::id(reads.BC)), BC=as.character(sread(reads.BC)))

##"CATTACGCGCTCGCGTAAGC" %in% names(frag.ranges.matched)
#' Extraction of fragments
#' ============================
#+ Extracting fragments.......


out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("-Xmx12g overwrite=true k=10 mink=18 rcomp=f qhdist=0 maskmiddle=t hammingdistance=1 findbestmatch=t ordered=t minlength=35 maxlength=75 threads=", detectCores(),
                      " in=", in.name.P7,
                      " out=", out.name.P7,
                      " lliteral=", "AGCAACCTCCAGAGAGGCAAC",
                      " rliteral=", "CAGACAAGCAGCTACCGCAGATGTCAACACACAAGGCGTTCTTCCAGGCATGGTCTGG", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t ,

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) # 

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")



in.name.P7 <- out.name.P7


source("retrieveFASTAQID.R")

FastQ1 <- readFastq(out.name.P5)
FastQ2 <- readFastq(out.name.P7)
FastQ1ID <- retrieveFASTAQID(FastQ1, PE=TRUE)
FastQ2ID <- retrieveFASTAQID(FastQ2, PE=TRUE)


hits <- intersect(FastQ2ID,FastQ1ID)

FastQ1Subset <- FastQ1[match(hits,FastQ1ID)]
FastQ2Subset <- FastQ2[match(hits,FastQ2ID)]

system(paste("mv ", out.name.P7, " ./data/fragments_", name.out, ".fastq.gz", sep=""))
system(paste("mv ", out.name.P5, " ./data/barcodes_", name.out, ".fastq.gz", sep=""))

unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) #Cleanup of temp files

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()