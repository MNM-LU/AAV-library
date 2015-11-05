require(ShortRead)
setwd("~/Dropbox (Bjorklund Lab)/Shared/NGS data")
source("R analysis/Functions/retrieveFASTAQID.R")

FastQFile1 <- "Original sequencing files/141114_M01551_16642633_000000000-AAEGT/010minusSalI1_S3_L001_R1_001.fastq.gz"
FastQFile2 <- "Original sequencing files/141114_M01551_16642633_000000000-AAEGT/010minusSalI2_S4_L001_R1_001.fastq.gz"
fileName1 <- sub(".fastq.gz", "", basename(FastQFile1))
fileName2 <- sub(".fastq.gz", "", basename(FastQFile2))

FastQ1 <- readFastq(FastQFile1)
FastQ2 <- readFastq(FastQFile2)
FastQ1ID <- retrieveFASTAQID(FastQ1, PE=TRUE)
FastQ2ID <- retrieveFASTAQID(FastQ2, PE=TRUE)


hits <- intersect(FastQ2ID,FastQ1ID)

FastQ1Subset <- FastQ1[match(hits,FastQ1ID)]
FastQ2Subset <- FastQ2[match(hits,FastQ2ID)]

writeFastq(FastQ1Subset,paste(fileName1, "_matched.fastq.gz", sep = ""), compress = TRUE)
writeFastq(FastQ2Subset,paste(fileName2, "_matched.fastq.gz", sep = ""), compress = TRUE)