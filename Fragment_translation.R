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


reads.trim <- readFastq("data/fragments_2015-11-04_AAVlibrary.fastq.gz")
reads.BC <- readFastq("data/barcodes_2015-11-04_AAVlibrary.fastq.gz")


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

reads.BC <- reads.BC[width(sread(reads.trim)) < 75L & width(sread(reads.trim)) > 38L]
reads.trim <- reads.trim[width(sread(reads.trim)) < 80L & width(sread(reads.trim)) > 35L]

source("retrieveFASTAQID.R")
FastQ1ID <- retrieveFASTAQID(reads.BC, PE=TRUE)
FastQ2ID <- retrieveFASTAQID(reads.trim, PE=TRUE)

hits <- intersect(FastQ2ID,FastQ1ID)

reads.BC <- reads.BC[match(hits,FastQ1ID)]
reads.trim <- reads.trim[match(hits,FastQ2ID)]
reads.table <- sread(reads.trim)
names(reads.table) <- sread(reads.BC)
reads.table.list <- split(reads.table, names(reads.table))
rm(reads.BC,reads.trim,FastQ1ID,FastQ2ID,hits,reads.table)
reads.table.list <- reads.table.list[1:1000]

reads.table.list <- reads.table.list[mcmapply(length,reads.table.list,mc.preschedule = TRUE, mc.cores = detectCores()) > 3]

reads.table.list <- reads.table.list[sapply(reads.table.list,length) > 3]

reads <- reads.table.list[[99]]


uniques <- rev(sort(table(reads)))
read.match <- adist(names(uniques),LUT.dna$Sequence)
names(uniques) <- unlist(lapply(1:nrow(read.match), function(i) which.min(read.match[i,])))
uniques <- split(uniques, names(uniques))
uniques <- unlist(lapply(uniques,sum))
uniques <- uniques[which.max(uniques)]
rowMins(read.match)
read.match <- sweep(read.match,MARGIN=1,uniques,`/`)


pmin(read.match)
read.match <- colMedians(read.match)
which.min(read.match)
lv.score <- min(read.match)/length(uniques)
matched.ref <- LUT.dna$Sequence[which.min(read.match)]
return(cbind(matched.ref,lv.score))

reads <- reads.table.list[[99]]
reads <- reads[4]
uniques <- rev(sort(table(reads)))
read.match <- adist(names(uniques),LUT.dna$Sequence)
LUT.dna$Sequence[which.min(read.match)]
which.min(read.match)




tmp %/% t(as.matrix(uniques))
read.match <- lapply(uniques,function(x) adist(names(x),LUT.dna$Sequence))

tmp <- read.match[[1]]
min(tmp/10)
match.pair <- function(x){
  matches <- adist(x,LUT.dna$Sequence)
  lv.score <- min(matches)
  matched.ref <- LUT.dna$Sequence[which.min(matches)]
  return(cbind(matched.ref,lv.score))
}
strt4<-Sys.time()
match.out <- mclapply(reads.sub, match.pair.c, mc.preschedule = TRUE, mc.cores = 32)
table.sub <- as.data.frame(do.call(rbind,match.out))

print("Total fragment translation time:")
print(Sys.time()-strt4)




print("Reduce pure.......")
reads.table.list <- mcmapply(reduce,reads.table.list, with.revmap=TRUE, mc.preschedule = TRUE, mc.cores = 32)
print(Sys.time()-strt2)



frag.ranges.subset.list <- split(frag.ranges.subset, names(frag.ranges.subset))



table.frag <- as.data.frame((rev(sort(table(sread(reads.trim)))))[1:100])
colnames(table.frag) <- c("Fragment and readcount")
knitr::kable(table.frag, format = "markdown")



reads.unique <- unique(sread(reads.trim))

## Find exact matches
match.matrix <- match(reads.unique, LUT.dna$Sequence)
reads.unique.match <- reads.unique[!is.na(match.matrix)]
reads.unique <- reads.unique[is.na(match.matrix)]

reads.sub <- reads.unique[1:100]
# strt4<-Sys.time()
# match.out1 <- lapply(reads.sub, function(x) which.min(stringdist(x,LUT.dna$Sequence, nthread = 32)))
# print(Sys.time()-strt4)


# match.out <- amatch(reads.sub, LUT.dna$Sequence, maxDist=6, nthread = 32)
# strt4<-Sys.time()
# match.out <- mclapply(reads.sub, function(x) LUT.dna$Sequence[which.min(adist(x,LUT.dna$Sequence))], mc.preschedule = TRUE, mc.cores = 32)
# reads.sub <- unlist(match.out)
# 
# print("Total fragment translation time:")
# print(Sys.time()-strt4)
# 
# match.out <- amatch(reads.sub, LUT.dna$Sequence, maxDist=6, nthread = 32)

match.pair <- function(x){
  matches <- adist(x,LUT.dna$Sequence)
  lv.score <- min(matches)
  matched.ref <- LUT.dna$Sequence[which.min(matches)]
  return(cbind(matched.ref,lv.score))
}
strt4<-Sys.time()
match.out <- mclapply(reads.sub, match.pair.c, mc.preschedule = TRUE, mc.cores = 32)
table.sub <- as.data.frame(do.call(rbind,match.out))

print("Total fragment translation time:")
print(Sys.time()-strt4)

table.match <- as.data.frame(reads.unique.match)
colnames(table.match) <- "V1"
table.match$V2 <- 0


reads.unique <- rbind(table.match,table.sub)

LUT.dna$Sequence[as.integer(names(reads.sub))]




reads.BC <- reads.BC[!is.na(match(sread(reads.trim), reads.unique$V1))]
reads.trim <- reads.trim[!is.na(match(sread(reads.trim), reads.unique$V1))]

sread(reads.trim) <- LUT.dna$Sequence[as.integer(names(reads.unique[match(sread(reads.trim), reads.unique)]))]



ShortReads::id(reads.BC) <- names(reads.sub[match(sread(reads.trim),reads.sub)])
names(frag.ranges.subset) <- barcodeTable$BC[match(names(frag.ranges.subset),barcodeTable$ID)]
table.BC.sc <- reads.sub[reads.unique %in% sread(reads.trim)]

#' Align fragments to reference
#' ============================
#+ Align to reference...

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --threads ",detectCores()," --local -D 20 -R 3 -N 1 -L 10 ",
                         "-i S,1,0.10 --ma 3 --no-unal --phred33 -x ",fragmentTemplate," -U ",out.name.P7, " -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE)

#   system(paste("samtools view -@ ",detectCores()," -Sh ", name.bowtie, ".sam | grep -v \"XS:i:\" > ",
#                  name.bowtie, "_filtered.sam",  sep = ""))
system(paste("samtools view -@ ",detectCores()," -bS ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

frag.ranges.ok <- frag.ranges[width(frag.ranges) > 38 & width(frag.ranges) < 75]
unique(granges(frag.ranges))

frag.ranges.subset <- frag.ranges[names(frag.ranges) %in% barcodeTable$ID]
names(frag.ranges.subset) <- barcodeTable$BC[match(names(frag.ranges.subset),barcodeTable$ID)]
frag.ranges.subset[names(frag.ranges.subset) %in% "GAGCATGGAAGCATGGCTGT"]

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