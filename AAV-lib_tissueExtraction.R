
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


#+ setup, include=FALSE
opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(comment = NA)


config <- read.table("config_tissue.txt", header = FALSE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
colnames(config) <- c("Parameter", "Value")
#setwd("~/Dropbox (Bjorklund Lab)/Shared/NGS data")
script.dir <- normalizePath("../SharedFunctions/")
# source(file.path(script.dir, "buildFragLengthTable.R"))
# source(file.path(script.dir, "makePEfastq.R"))
# source(file.path(script.dir, "retrieveFASTAQID.R"))
source(file.path(script.dir, "buildFragTableSC3.R"))

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

load("singleContfragments.rda")
output.Table.single <- output.Table
load("multipleContfragments.rda")
output.Table <- rbind(output.Table, output.Table.single)
load("alignedLibraries.rda")
load("LUTdna.rda")

output.Table <- na.omit(output.Table)
nrow(output.Table)
#output.Table <-output.Table[(output.Table$LV < 4),] #(output.Table$mCount/output.Table$tCount > 0.5) & 
nrow(output.Table)






#' Selection of real amplicons
#' ============================
#+ Selecting real amplicons.......

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("-Xmx12g overwrite=true k=10 rcomp=f skipr1=t qhdist=0 maskmiddle=t hammingdistance=0 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "CGCCACAACATCGAGGACGGCAGCGTG", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t , ATATCATGGCCGACAAGCAGA

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)


colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7

nr.Reads <- as.integer(system(paste("gunzip -c ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE)) #Stores the read count utilized
print(paste("Utilized sequences:", nr.Reads))


#' Extraction of barcodes
#' ============================
#+ Extracting barcodes.......

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
command.args <- paste("-Xmx12g overwrite=true k=12 rcomp=f skipr2=t qhdist=0 maskmiddle=t hammingdistance=2 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "ATAACTTCGTATA", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t ,

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)


colnames(sys.out) <- c("bbduk2 Identification of real barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7


out.name.BC <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=12 mink=12 hammingdistance=2 findbestmatch=t ",
                        "trd=t rcomp=f skipr2=t findbestmatch=f qhdist=0 minavgquality=0 ordered=t maxns=0 minlength=18 ",
                        "maxlength=22 threads=", detectCores()," in=", shQuote(in.name.P5),
                        " out=", out.name.BC,
                        " lliteral=", "GGCCTAGCGGCCGCTTTACTT",
                        " rliteral=", "ATAACTTCGTATA",
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) #" fliteral=",id.uncut,
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")
rm(sys.out)

reads.BC <- readFastq(out.name.BC)
sread(reads.BC)
unique(sread(reads.BC))
barcodeTable <- data.table(ID=as.character(ShortRead::id(reads.BC)), BC=as.character(sread(reads.BC)))

setkey(barcodeTable,BC)

#' Starcode based barcode reduction
#' ============================
#+ Reducing barcodes.......

out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")
system(paste("gunzip -c ",out.name.BC," | starcode -t ",detectCores()/2," --print-clusters -d",
             config$Value[8]," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
       intern = TRUE, ignore.stdout = FALSE)

table.BC.sc <- read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                          stringsAsFactors = FALSE, fill=FALSE)
table.BC.sc$V2 <- NULL
list.BC.sc <- split(table.BC.sc, rownames(table.BC.sc))
list.BC.sc <- lapply(list.BC.sc, function(x) strsplit(as.character(x), ","))
list.BC.sc <- lapply(list.BC.sc, rbind)
table.BC.sc <- data.table(cbind(unlist(list.BC.sc, use.names = TRUE)), keep.rownames=TRUE)
invisible(table.BC.sc[,rn:=gsub("[0-9]","",rn)])
droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
print(paste("Dropped BCs in Starcode:", droppedBC))

setnames(table.BC.sc,'V1','BC')
setnames(table.BC.sc,'rn','scBC')

setkey(table.BC.sc,BC)
table.BC.sc<- unique(table.BC.sc)
BCs.sc.counts <- table(table.BC.sc$BC)
BCs.sc.counts.single <- BCs.sc.counts[BCs.sc.counts == 1]
table.BC.sc <- table.BC.sc[table.BC.sc$BC %in% names(BCs.sc.counts.single)]

barcodeTable <- merge(barcodeTable,table.BC.sc, by="BC", all = FALSE, all.x = FALSE)
rm(table.BC.sc)
rm(reads.BC)
setnames(barcodeTable,'BC','oldBC')
setnames(barcodeTable,'scBC','BC')
setkey(barcodeTable,BC)

RetainedBC <- length(unique(barcodeTable$oldBC))
scBC <- length(unique(barcodeTable$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))


table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$oldBC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
knitr::kable(table.frag, format = "markdown")

table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$BC))))[1:15]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
knitr::kable(table.frag, format = "markdown")

invisible(barcodeTable[,oldBC:=NULL])


BCcount <- rev(sort(table(barcodeTable$BC)))

foundFrags <- output.Table[match(names(BCcount), output.Table$BC),]
foundFrags$RNAcount <- as.integer(BCcount)
foundFrags <- na.omit(foundFrags)
foundFrags$fragment <- LUT.dna$Sequence[as.integer(foundFrags$LUTnr)]

matchRange <- function(idxFrag) {
  #idxFrag <- 23
  machRanges <- which(names(allFragments.ranges) == foundFrags$fragment[idxFrag])
  return(cbind(machRanges,idxFrag))
}
match.ranges.list <- mclapply(1:nrow(foundFrags), matchRange, mc.preschedule = TRUE, mc.cores = detectCores())
match.ranges <- do.call(rbind, match.ranges.list)
foundFragments.ranges <- allFragments.ranges[match.ranges[,1]]
mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges), foundFrags[match.ranges[,2],2:6])

o = order(-mcols(foundFragments.ranges)$RNAcount)
foundFragments.ranges <- foundFragments.ranges[o]
foundFragments.ranges[1:10]

assign(paste("found.",name.out, sep=""), foundFragments.ranges)
save(list = paste("found.",name.out, sep=""), file=paste("output/","found.",name.out,".rda", sep=""))
unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) #Cleanup of temp files

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()