
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


config <- read.table("config.txt", header = FALSE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
colnames(config) <- c("Parameter", "Value")
#setwd("~/Dropbox (Bjorklund Lab)/Shared/NGS data")
script.dir <- "~/Dropbox (Bjorklund Lab)/Shared/NGS data/R analysis/Functions/"
source(file.path(script.dir, "buildFragLengthTable.R"))
source(file.path(script.dir, "makePEfastq.R"))
source(file.path(script.dir, "retrieveFASTAQID.R"))
source(file.path(script.dir, "buildFragTableSC3.R"))

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
output.table$Reads <- as.integer(system(paste("zcat ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE)) #Stores the read count utilized
print(paste("Utilized sequences:", output.table$Reads[1]))


#' Extraction of barcodes
#' ============================
#+ Extracting barcodes.......
out.name.BC <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=12 hammingdistance=1 findbestmatch=t ",
                        "trd=t rcomp=f findbestmatch=f qhdist=0 minavgquality=30 maxns=0 minlength=18 ",
                        "maxlength=22 threads=", detectCores()," in=", shQuote(in.name.P5), 
                        " out=", out.name.BC," lref=", id.BC.L,
                        " rref=", id.BC.R," fliteral=",id.uncut,
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE)

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")
rm(sys.out)

reads.BC <- readFastq(out.name.BC)
sread(reads.BC)
output.table$OrigBC <- length(unique(sread(reads.BC)))
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
output.table$droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1) %in% unique(sread(reads.BC)))
print(paste("Dropped BCs in Starcode:", output.table$droppedBC[1]))

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

output.table$RetainedBC <- length(unique(barcodeTable$oldBC))
output.table$scBC <- length(unique(barcodeTable$BC))
print(paste("Original unique barcodes:", output.table$RetainedBC[1]))
print(paste("SC reduced unique barcodes:", output.table$scBC[1]))


table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$oldBC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("Original BC", "Count"))
knitr::kable(table.frag, format = "markdown")

table.frag <- data.table(as.data.frame((rev(sort(table(barcodeTable$BC))))[1:10]), keep.rownames=TRUE)
setnames(table.frag, colnames(table.frag), c("SC reduced BC", "Count"))
knitr::kable(table.frag, format = "markdown")

invisible(barcodeTable[,oldBC:=NULL])

##"CATTACGCGCTCGCGTAAGC" %in% names(frag.ranges.matched)
#' Extraction of fragments
#' ============================
#+ Extracting fragments.......


# reads.untrim <- readFastq(in.name.P5)
# reads.untrim <- reads.untrim[width(sread(reads.untrim))>= 90]
# 
# pattern.pre <- "GTATGTTGTTCTGG"
# found.match <- vcountPattern(pattern.pre, sread(reads.untrim),
#                              max.mismatch=0, min.mismatch=0,
#                              with.indels=FALSE, fixed=TRUE,
#                              algorithm="auto")
# reads.with3p <- reads.untrim[as.logical(found.match)]
# rev(sort(table(sread(reads.with3p))))[1:10]
# sum(found.match)
# 
# 
# 
# 
# pattern.pre <- "CATGGACGAGCTGTACAAGT"
# found.match <- vcountPattern(pattern.pre, sread(reads.untrim),
#                              max.mismatch=0, min.mismatch=0,
#                              with.indels=FALSE, fixed=TRUE,
#                              algorithm="auto")
# reads.with3p <- reads.untrim[as.logical(found.match)]
# sum(found.match)
# rev(sort(table(sread(reads.with3p))))[1:10]
# 
# 
# unique(sread(reads.with3p))
# 
# pattern.post <- "CAGACAAGCAGCTACCGCA"
# found.match <- vcountPattern(pattern.post, sread(reads.untrim),
#               max.mismatch=1, min.mismatch=0,
#               with.indels=TRUE, fixed=TRUE,
#               algorithm="auto")
# 
# sum(found.match)


out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=18 mink=10 qhdist=0 minlength=44 maxlength=75 hammingdistance=1 findbestmatch=t threads=",detectCores(),
                        " in=", in.name.P7, 
                        " out=", out.name.P7,
                        " lref=", id.backbone.L,
                        " rref=", id.backbone.R, 
                        " fliteral=",id.uncut,
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t 
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

reads.trim <- readFastq(out.name.P7)

table.frag <- as.data.frame((rev(sort(table(sread(reads.trim)))))[1:100])
colnames(table.frag) <- c("Fragment and readcount")
knitr::kable(table.frag, format = "markdown")


#' Align fragments to reference
#' ============================
#+ Align to reference...

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")
if (paired.alignment){
  sys.out <-  system(paste("bowtie2 --no-discordant --no-mixed --fr --threads ",detectCores()
                           ,"  --local --very-sensitive-local --no-unal --gbar 10 --rdg 10,20 ",
                           "--minins 20 --maxins 2500 --phred33 -x ",
                           fragmentTemplate," -1 ",out.name.P5, " -2 ", out.name.P7, " -S ", 
                           name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE)
  
  system(paste("samtools view -@ ",detectCores()," -Shf 0x2 ", name.bowtie, ".sam | grep -v \"XS:i:\" > ",
               name.bowtie, "_filtered.sam",  sep = ""))
  system(paste("samtools view -@ ",detectCores()," -bS ", name.bowtie, "_filtered.sam > ",name.bowtie, 
               ".bam",  sep = ""))
  system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",name.bowtie, "_sort",  sep = ""))
  
  frag.ranges <- readGAlignmentPairs(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
  
} else {
  if (align.p7) {
    sys.out <-  system(paste("bowtie2 --threads ",detectCores()," --local -D 20 -R 3 -N 1 -L 10 ",
                             "-i S,1,0.10 --ma 3 --no-unal --phred33 -x ",fragmentTemplate," -U ",out.name.P7, " -S ", 
                             name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE)
  } else{
    sys.out <-  system(paste("bowtie2 --threads ",detectCores()," --local -D 20 -R 3 -N 1 -L 10 ",
                             "-i S,1,0.10 --ma 3 --no-unal --phred33 -x ",fragmentTemplate," -U ",out.name.P5, " -S ", 
                             name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE)
  }
  
  system(paste("samtools view -@ ",detectCores()," -bS ", name.bowtie, ".sam > ",
               name.bowtie, ".bam",  sep = ""))
  system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
               name.bowtie, "_sort",  sep = ""))
  
  frag.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
}
sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Extraction of barcodes")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

frag.ranges <- granges(frag.ranges)
frag.ranges <- frag.ranges[seqnames(frag.ranges) != "backbone"]


#' Generation of summary table
#' ============================
#+ Generating summary table.......

frag.ranges.matched <- buildFragTable(frag.ranges, barcodeTable, run.subset=FALSE, max.cores)

save(frag.ranges.matched, file = paste(name.out, "_fragRanges.rda", sep = ""), compress=TRUE)

frag.ranges.matched <- frag.ranges.matched[(!(frag.ranges.matched$count >1 & 
                                                (frag.ranges.matched$hitsStart==1 | frag.ranges.matched$hitsEnd==1)))]

num.size.seq  <- max(max(width(frag.ranges)), as.numeric(frag.ranges@seqinfo@seqlengths))

#' Plot barcode count
#' ==================
#+ GeneratingPlots
uniqueBCs <- data.frame(BCcount = frag.ranges.matched$count,row.names=names(frag.ranges.matched))
num.count.max.BC <- max(uniqueBCs)
num.size.bin  <- ceiling(num.count.max.BC/100)
p.linear <- ggplot(uniqueBCs, aes(uniqueBCs$BCcount)) + 
  geom_histogram(binwidth = num.size.bin, fill=NA, color="black") + theme_bw() + 
  labs(title = "Barcode count distribution", y = "Number of barcodes", x = "reads per barcode") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4)) +
  scale_x_continuous(limit=c(0,num.count.max.BC), expand =c(0,0))

p.log <-ggplot(uniqueBCs, aes(uniqueBCs$BCcount)) + 
  geom_histogram(binwidth = num.size.bin, fill=NA, color="black") + theme_bw() + 
  scale_y_continuous(trans=log2_trans()) + 
  labs(title = "Barcode count distribution", y = "Number of barcodes (log2)", x = "reads per barcode") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4)) +
  scale_x_continuous(limit=c(0,num.count.max.BC), expand =c(0,0))

plot.montage <- tracks("Linear BC count" = p.linear, "log2 BC count" = p.log, heights = c(8, 8), 
                       title="Barcode count distribution") + ylab("") + theme_tracks_sunset()

suppressWarnings(print(plot.montage)) #Removes a warning that 0 is not defined after log-transformed


opts_chunk$set(fig.width = 5.5, fig.height = 6)
#' Make a beanplot
#' ===============
#+ beanplot
frag.ranges.multi  <- frag.ranges.matched[frag.ranges.matched$count >1,]

output.table$meanCount <- mean(frag.ranges.multi$hits/frag.ranges.multi$count)
output.table$sdCount <- sd(frag.ranges.multi$hits/frag.ranges.multi$count)
output.table$meanWidth <- mean(width(frag.ranges.multi))
output.table$sdWidth <- sd(width(frag.ranges.multi))
output.table[c("q0","q25","q50","q75","q100")] <- quantile(width(frag.ranges.multi))
output.table$Name <- as.character(seqnames(frag.ranges.multi)[1])

print(paste("Valid barcodes ulilized:", length(frag.ranges.multi)))
frag.ranges.multi.fwd  <- frag.ranges.multi[strand(frag.ranges.multi) == "+"]
frag.ranges.multi.rev  <- frag.ranges.multi[strand(frag.ranges.multi) == "-"]
purity.P5 <- as.data.frame(append((frag.ranges.multi.fwd$hitsStart/frag.ranges.multi.fwd$count),
                                  (frag.ranges.multi.rev$hitsEnd/frag.ranges.multi.rev$count)))
purity.P7 <- as.data.frame(append((frag.ranges.multi.fwd$hitsEnd/frag.ranges.multi.fwd$count),
                                  (frag.ranges.multi.rev$hitsStart/frag.ranges.multi.rev$count)))

output.table$meanP5 <- mean(purity.P5[[1]])
output.table$sdP5 <- sd(purity.P5[[1]])
output.table$meanP7 <- mean(purity.P7[[1]])
output.table$sdP7 <- sd(purity.P7[[1]])

purity.P5[,2]  <- paste(name.out, " 1", sep = "")
purity.P7[,2]  <- paste(name.out, " 2", sep = "")
colnames(purity.P5)[1:2]  <- c("pure.fraction", "Sample")
colnames(purity.P7)[1:2]  <- c("pure.fraction", "Sample")
purity  <- rbind(purity.P5,purity.P7)

beanplot(pure.fraction ~ Sample, data = purity, ll = 0.04, what=c(1,1,1,0), bw = "nrd0", log="", 
         main = "Recombination analysis", ylab = "Recombination fraction [1 = no recombination]", 
         side = "both", border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black", "grey"), legend = c("P5 end", "P7 end"))

purity.total <- as.data.frame(frag.ranges.multi$hits/frag.ranges.multi$count)
purity.total[,2]  <- paste(name.out, " 1", sep = "")
colnames(purity.total)[1:2]  <- c("pure.fraction", "Sample")

beanplot(pure.fraction ~ Sample, data = purity.total, ll = 0.04, what=c(1,1,1,0), bw = "nrd0", log="", 
         main = "Total recombination analysis", ylab = "Recombination fraction [1 = no recombination]", 
         border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black"), legend = c("Total"))

purity.total <- as.data.frame(rowMax(cbind(frag.ranges.multi$hitsStart, frag.ranges.multi$hitsEnd))/frag.ranges.multi$count)
purity.total[,2]  <- paste(name.out, " 1", sep = "")
colnames(purity.total)[1:2]  <- c("pure.fraction", "Sample")

beanplot(pure.fraction ~ Sample, data = purity.total, ll = 0.04, what=c(1,1,1,0), bw = "nrd0", log="", 
         main = "Best end recombination analysis", ylab = "Recombination fraction [1 = no recombination]", 
         border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black"), legend = c("Best end"))

frag.ranges.unique  <- unique(frag.ranges.multi)

opts_chunk$set(fig.width = 7.5, fig.height = 8)
#'Plot coverage
#'=================
p.all <- suppressMessages(autoplot(frag.ranges, stat="coverage"))
p.matched <- suppressMessages(autoplot(frag.ranges.matched, stat="coverage"))
p.unique <- suppressMessages(autoplot(frag.ranges.unique, stat="coverage"))

plot.montage <- tracks("All fragments" = p.all, "Fragments w. unique BCs" = p.matched, 
                       "Unique fragments" = p.unique, heights = c(5, 5, 5), title="Coverage plots") + 
  ylab("") + theme_tracks_sunset()
print(plot.montage)

#'Plot fragment distribution
#'=================
#+ fragmentDistribution

p.all <- ggplot(frag.ranges, aes(width(frag.ranges))) + 
  geom_histogram(binwidth = 10, fill=NA, color="black") +
  scale_x_continuous(limit=c(0,num.size.seq), expand =c(0,0))

p.matched <- ggplot(frag.ranges.matched, aes(width(frag.ranges.matched))) + 
  geom_histogram(binwidth = 10, fill=NA, color="black") +
  scale_x_continuous(limit=c(0,num.size.seq), expand =c(0,0))

p.unique <- ggplot(frag.ranges.unique, aes(width(frag.ranges.unique)), xmin = 0, xmax=num.size.seq) + 
  geom_histogram(binwidth = 10, fill=NA, color="black") +
  scale_x_continuous(limit=c(0,num.size.seq), expand =c(0,0))

plot.montage <- tracks("All fragments" = p.all, "Fragments w. unique BCs" = p.matched, 
                       "Unique fragments" = p.unique, heights = c(8, 8, 8), title="Fragment length distribution") + 
  ylab("") + theme_tracks_sunset()
print(plot.montage)

#' Generate unique fragment plots
#' ==============================
#+ uniqueFragments
o = order(start(frag.ranges.unique),end(frag.ranges.unique))
frag.ranges.unique <- frag.ranges.unique[o]
frag.ranges.uniqueFwd  <- frag.ranges.unique[strand(frag.ranges.unique) == "+"]
frag.ranges.uniqueRev  <- frag.ranges.unique[strand(frag.ranges.unique) == "-"]
output.table$uniqueFragFwd <- length(frag.ranges.uniqueFwd)
output.table$uniqueFragRev <- length(frag.ranges.uniqueRev)
print(paste("Unique Fragments FWD:", output.table$uniqueFragFwd[1]," REV:", output.table$uniqueFragRev[1]))

p.fwd <- ggplot(frag.ranges.uniqueFwd) + ggplot2::geom_rect(aes(xmin = start, ymin = 1:length(frag.ranges.uniqueFwd), 
                                                                xmax = end, ymax = 1:length(frag.ranges.uniqueFwd) + 1))

p.rev <- ggplot(frag.ranges.uniqueRev) + ggplot2::geom_rect(aes(xmin = start, ymin = 1:length(frag.ranges.uniqueRev), 
                                                                xmax = end, ymax = 1:length(frag.ranges.uniqueRev) + 1))
plot.montage <- tracks("+ strand" = p.fwd, "- strand" = p.rev, heights = c(8, 8), 
                       title="Fragment plots", xlab="Base of reference sequence") + 
  ylab("") + theme_tracks_sunset()
print(plot.montage)

p.count.matched <- ggplot(frag.ranges.matched) + ggplot2::geom_rect(aes(xmin = start, ymin = count, xmax = end, ymax = count + 1))
p.count.unique <- ggplot(frag.ranges.unique) + ggplot2::geom_rect(aes(xmin = start, ymin = count, xmax = end, ymax = count + 1))
plot.montage <- tracks("All matched" = p.count.matched, "Unique fragments" = p.count.unique, heights = c(8, 8), 
                       title="Fragment counts", xlab="Base of reference sequence") + 
  ylab("Read count") + theme_tracks_sunset()
print(plot.montage)


knitr::kable(output.table[,1:8], format = "markdown")
knitr::kable(output.table[,9:13], format = "markdown")
knitr::kable(output.table[,14:as.integer(ncol(output.table))], format = "markdown")

write.csv(output.table, file = paste(name.out, "_outTable.csv", sep = ""))

unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) #Cleanup of temp files

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()