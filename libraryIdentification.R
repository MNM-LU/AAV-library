suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

LUT.dna <- read.table("Complete fragment list for Custom array 2015-02-10.txt", header = TRUE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
LUT.dna <- data.table(LUT.dna)
invisible(LUT.dna[,Sequence:=gsub("aacctccagagaggcaac","",Sequence)])
invisible(LUT.dna[,Sequence:=gsub("agacaagcagctaccgca","",Sequence)])
invisible(LUT.dna[,Sequence:=toupper(Sequence)])
setkey(LUT.dna, "Sequence")
LUT.dna <- unique(LUT.dna)
LUT.dna$Names <- LUT.dna$Sequence

#output.Table$LUTseq <- LUT.dna$Sequence[as.numeric(output.Table$LUTnr)]
LUT.14aaG4S <- LUT.dna[substr(LUT.dna$Sequence,1,15) == "GGAGGCGGAGGAAGT"]
LUT.remaining <- LUT.dna[!(substr(LUT.dna$Sequence,1,15) == "GGAGGCGGAGGAAGT")]
LUT.14aaA5 <- LUT.remaining[substr(LUT.remaining$Sequence,1,15) == "GCTGCTGCAGCAGCC"]
LUT.remaining <- LUT.remaining[!(substr(LUT.remaining$Sequence,1,15) == "GCTGCTGCAGCAGCC")]
LUT.22aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 72L & substr(LUT.remaining$Sequence,1,3) == "GCT"]
LUT.remaining <- LUT.remaining[!(nchar(LUT.remaining$Sequence) == 72L & substr(LUT.remaining$Sequence,1,3) == "GCT")]
LUT.14aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 48L & substr(LUT.remaining$Sequence,1,3) == "GCT"]

LUT.14aa$Sequence <- substr(LUT.14aa$Sequence,4,45)
LUT.14aaG4S$Sequence <- substr(LUT.14aaG4S$Sequence,16,57)
LUT.14aaA5$Sequence <- substr(LUT.14aaA5$Sequence,16,57)
LUT.22aa$Sequence <- substr(LUT.22aa$Sequence,4,69)

LUT.14aa.seq = ShortRead(DNAStringSet(LUT.14aa$Sequence), BStringSet(LUT.14aa$Names))
writeFasta(LUT.14aa.seq,"LUT.14aa.fa")

LUT.14aaG4S.seq = ShortRead(DNAStringSet(LUT.14aaG4S$Sequence), BStringSet(LUT.14aaG4S$Names))
writeFasta(LUT.14aaG4S.seq,"LUT.14aaG4S.fa")

LUT.14aaA5.seq = ShortRead(DNAStringSet(LUT.14aaA5$Sequence), BStringSet(LUT.14aaA5$Names))
writeFasta(LUT.14aaA5.seq,"LUT.14aaA5.fa")

LUT.22aa.seq = ShortRead(DNAStringSet(LUT.22aa$Sequence), BStringSet(LUT.22aa$Names))
writeFasta(LUT.22aa.seq,"LUT.22aa.fa")



#' Align fragments to reference
#' ============================
#+ Align to reference...

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores()," --local --score-min 'C,0,-1' -f -a",
                         " -x bowtieIdx/libIdx -U LUT.14aa.fa -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE) #  --ma 3 --no-unal --phred33  -D 20 -R 3 -N 1 -L 10 -i S,1,0.10


sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -b ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aa.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aa.ranges))
length(unique(names(frag14aa.ranges)))
length(unique(LUT.14aa$Sequence))

#14aaG4S

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores()," --local --score-min 'C,0,-1' -f -a",
                         " -x bowtieIdx/libIdx -U LUT.14aaG4S.fa -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE) #  --ma 3 --no-unal --phred33  -D 20 -R 3 -N 1 -L 10 -i S,1,0.10


sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -b ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aaG4S.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aaG4S.ranges))
length(unique(names(frag14aaG4S.ranges)))
length(unique(LUT.14aaG4S$Sequence))

#14aaA5

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores()," --local --score-min 'C,0,-1' -f -a",
                         " -x bowtieIdx/libIdx -U LUT.14aaA5.fa -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE) #  --ma 3 --no-unal --phred33  -D 20 -R 3 -N 1 -L 10 -i S,1,0.10


sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -b ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag14aaA5.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag14aaA5.ranges))
length(unique(names(frag14aaA5.ranges)))
length(unique(LUT.14aaA5$Sequence))

#22aa

name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2 --non-deterministic --threads ",detectCores()," --local --score-min 'C,0,-1' -f -a",
                         " -x bowtieIdx/libIdx -U LUT.22aa.fa -S ", 
                         name.bowtie, ".sam 2>&1",  sep = ""), intern = TRUE, ignore.stdout = FALSE) #  --ma 3 --no-unal --phred33  -D 20 -R 3 -N 1 -L 10 -i S,1,0.10


sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("Bowtie 2 alignment to library")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[1:lengthOut,], format = "markdown")

system(paste("samtools view -@ ",detectCores()," -b ", name.bowtie, ".sam > ",
             name.bowtie, ".bam",  sep = "")) 
system(paste("samtools sort -@ ",detectCores()," ", name.bowtie, ".bam ",
             name.bowtie, "_sort",  sep = ""))

frag22aa.ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names=TRUE)
length(names(frag22aa.ranges))
length(unique(names(frag22aa.ranges)))
length(unique(LUT.22aa$Sequence))
mcols(frag14aa.ranges)$structure <- "14aa"
mcols(frag22aa.ranges)$structure <- "22aa"
mcols(frag14aaA5.ranges)$structure <- "14aaA5"
mcols(frag14aaG4S.ranges)$structure <- "14aaG4S"
allFragments.ranges <- append(frag14aa.ranges,frag22aa.ranges)
allFragments.ranges <- append(allFragments.ranges,frag14aaA5.ranges)
allFragments.ranges <- append(allFragments.ranges,frag14aaG4S.ranges)

save(allFragments.ranges, file="alignedLibraries.rda")
