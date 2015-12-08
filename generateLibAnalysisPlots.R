suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(beanplot))

load("~/Dropbox (Bjorklund Lab)/R-projects/AAV-library/completeLibraryRanges.rda")
purity.table <- data.table(mcols(complete.ranges)$mCount/mcols(complete.ranges)$tCount)
purity.table$BCwidth <- width(mcols(complete.ranges)$BC)
qplot(mcols(complete.ranges)$tCount, stat = "ecdf", geom = "step")
beanplot(data = purity.table$V1, ll = 0.04, what=c(1,1,1,0), bw = "nrd0", log="", 
         main = "Best end recombination analysis", ylab = "Recombination fraction [1 = no recombination]", 
         border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black"), legend = c("Best end"))

qplot(mcols(complete.ranges)$tCount, stat = "ecdf", geom = "step") + 
  scale_x_continuous(limit=c(1,25), breaks=c(seq(1,50,10)), expand =c(0,0))

fill.values <- c("A" = rgb(193,210,234, maxColorValue = 255), "B" = rgb(157,190,217, maxColorValue = 255), "c" = rgb(38,64,135, maxColorValue = 255), "D" = rgb(157,190,217, maxColorValue = 255),"E" = rgb(193,210,234, maxColorValue = 255))
names(fill.values) <- c("18","19","20","21","22")

ggplot(purity.table,aes(x = BCwidth, fill = as.character(BCwidth)))+geom_histogram(binwidth=1, origin = -0.5)+theme_bw() +
  scale_fill_manual(name = "Width", values = fill.values) +
  scale_x_continuous(limit=c(17,23), breaks=c(seq(18,22,1)), expand =c(0,0))

#Plot Venn diagrams of fragments
load("~/Dropbox (Bjorklund Lab)/R-projects/AAV-library/completeLibraryRanges.rda")
LUT.dna <- read.table("Complete fragment list for Custom array 2015-02-10.txt", header = TRUE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
LUT.dna <- data.table(LUT.dna)
invisible(LUT.dna[,Sequence:=gsub("aacctccagagaggcaac","",Sequence)])
invisible(LUT.dna[,Sequence:=gsub("agacaagcagctaccgca","",Sequence)])
invisible(LUT.dna[,Sequence:=toupper(Sequence)])
setkey(LUT.dna, "Sequence")
LUT.dna <- unique(LUT.dna)
LUT.dna$Names <- LUT.dna$Sequence


total.library <- readRDS("output/total.infectiveLib.rds")
total.str <- GAlignmentsList(readRDS("output/found.RatNr15_1000x_Str-15_RatNr19_1000x_Str-22_RatNr20_1000x_Str-24_RatNr21_1000x_Str-19.rds"),readRDS("output/found.RatNr1_100x_Str-7_RatNr7_100x_Str-3_RatNr8_100x_Str-11.rds"))
total.str <- cbind(unlist(total.str))[[1]]
total.SNctx <- GAlignmentsList(readRDS("output/found.RatNr1_100x_Ctx-6_RatNr7_100x_Ctx-2_RatNr8_100x_Ctx-10.rds"),readRDS("output/found.RatNr1_100x_SN-5_RatNr7_100x_SN-1_RatNr8_100x_SN-9.rds"),readRDS("output/found.RatNr15_1000x_Ctx-14_RatNr19_1000x_Ctx-21_RatNr21_1000x_Ctx-18.rds"),readRDS("output/found.RatNr15_1000x_SN-13_RatNr21_1000x_SN-17.rds"))
total.SNctx <- cbind(unlist(total.SNctx))[[1]]



seq.arry <- LUT.dna$Sequence
seq.lib <- unique(names(complete.ranges))
seq.AAV <- unique(names(total.library))
seq.str <- unique(names(total.str))
seq.SNCtx <- unique(names(total.SNctx))

venn.area1 <- length(seq.arry)
venn.area2 <- length(seq.lib)
venn.area3 <- length(seq.AAV)
venn.area3 <- length(seq.str)
venn.area3 <- length(seq.SNCtx)


isect.Str_SN <- length(intersect(seq.str, seq.SNCtx))

venn.n12 <- length(intersect(seq.arry,seq.lib))
venn.n23 <- length(intersect(seq.lib,seq.AAV))
venn.n13 <- length(intersect(seq.arry,seq.AAV))
venn.n123 <- length(intersect(intersect(seq.arry,seq.lib),seq.AAV))


output.table <- data.frame(NameArray=character(),
                           NameLib=character(),
                           NameAAV=character(),
                           NameStr=character(),
                           NameSN=character(),
                           ArrayStart=numeric(), 
                           ArrayEnd=numeric(), 
                           LibStart=numeric(), 
                           LibEnd=numeric(), 
                           AAVStart=numeric(),
                           AAVend=numeric(),
                           StrStart=numeric(),
                           StrEnd=numeric(),
                           SNStart=numeric(),
                           SNEnd=numeric(),
                           stringsAsFactors=FALSE) 
output.table[1:3,1] <- c("Array","None","None")
output.table[1:3,2] <- c("Lib","None","None")
output.table[1:3,3] <- c("AAV","None","None")
output.table[1:3,4] <- c("Str","None","None")
output.table[1:3,5] <- c("None","SN","None")
output.table[1:3,6:15] <- 0
output.table$ArrayEnd[1] <- output.table$LibEnd[2] <- output.table$AAVend[2] <- output.table$StrEnd[2] <- output.table$SNEnd[3] <- length(seq.arry)
output.table$LibStart[2] <- output.table$LibEnd[1] <- length(intersect(seq.arry,seq.lib))
output.table$AAVStart[2] <- output.table$AAVend[1] <- length(intersect(seq.lib,seq.AAV))
output.table$StrStart[2] <- output.table$StrEnd[1] <- length(intersect(seq.AAV, seq.str))
output.table$SNStart[2] <- output.table$SNEnd[1] <- length(seq.str)-length(intersect(seq.str, seq.SNCtx))
output.table$SNStart[3] <- output.table$SNEnd[2] <- output.table$SNStart[2] + length(seq.SNCtx)


fill.values <- c("Array" = rgb(193,210,234, maxColorValue = 255), "Lib" = rgb(157,190,217, maxColorValue = 255), "AAV" = rgb(38,64,135, maxColorValue = 255), "Str" = rgb(157,190,217, maxColorValue = 255),"SN" = rgb(193,210,234, maxColorValue = 255), "None" = rgb(255,255,255, maxColorValue = 255, alpha = 0))


ggplot(output.table) + 
  scale_x_continuous(limit=c(0,10), breaks=c(seq(1,9,2)), expand =c(0,0)) + 
  scale_fill_manual(name = "Library", values = fill.values) +
  theme(aspect.ratio=1) + 
  geom_rect(data=output.table, aes(fill=NameSN, ymax=SNEnd, ymin=SNStart, xmax=10, xmin=8)) +
  geom_rect(data=output.table, aes(fill=NameStr, ymax=StrEnd, ymin=StrStart, xmax=8, xmin=6)) +
  geom_rect(data=output.table, aes(fill=NameAAV, ymax=AAVend, ymin=AAVStart, xmax=6, xmin=4)) +
  geom_rect(data=output.table, aes(fill=NameLib, ymax=LibEnd, ymin=LibStart, xmax=4, xmin=2)) +
  geom_rect(data=output.table, aes(fill=NameArray, ymax=ArrayEnd, ymin=ArrayStart, xmax=2, xmin=0)) +
  coord_polar(theta="y") 





total.100x <- GAlignmentsList(readRDS("output/found.Cells293Nr3_100x_cDNA-27.rds"),
                             readRDS("output/found.primNeuronsNr7_100x_cDNA-29.rds"),
                             readRDS("output/found.RatNr1_100x_Ctx-6.rds"),
                             readRDS("output/found.RatNr1_100x_SN-5.rds"),
                             readRDS("output/found.RatNr1_100x_Str-7.rds"),
                             readRDS("output/found.RatNr1_100x_Th-8.rds"),
                             readRDS("output/found.RatNr7_100x_Ctx-2.rds"),
                             readRDS("output/found.RatNr7_100x_SN-1.rds"),
                             readRDS("output/found.RatNr7_100x_Str-3.rds"),
                             readRDS("output/found.RatNr7_100x_Th-4.rds"),
                             readRDS("output/found.RatNr8_100x_Ctx-10.rds"),
                             readRDS("output/found.RatNr8_100x_SN-9.rds"),
                             readRDS("output/found.RatNr8_100x_Str-11.rds"),
                             readRDS("output/found.RatNr8_100x_Th-12.rds"))
total.100x <- cbind(unlist(total.100x))[[1]]

total.1000x <- GAlignmentsList(readRDS("output/found.Cells293Nr2_1000x_cDNA-26.rds"),
                              readRDS("output/found.primNeuronsNr6_1000x_cDNA-28.rds"),
                              readRDS("output/found.RatNr15_1000x_Ctx-14.rds"),
                              readRDS("output/found.RatNr15_1000x_SN-13.rds"),
                              readRDS("output/found.RatNr15_1000x_Str-15.rds"),
                              readRDS("output/found.RatNr15_1000x_Th-16.rds"),
                              readRDS("output/found.RatNr19_1000x_Ctx-21.rds"),
                              readRDS("output/found.RatNr19_1000x_Str-22.rds"),
                              readRDS("output/found.RatNr19_1000x_Th-23.rds"),
                              readRDS("output/found.RatNr20_1000x_Str-24.rds"),
                              readRDS("output/found.RatNr20_1000x_Th-25.rds"),
                              readRDS("output/found.RatNr21_1000x_Ctx-18.rds"),
                              readRDS("output/found.RatNr21_1000x_SN-17.rds"),
                              readRDS("output/found.RatNr21_1000x_Str-19.rds"),
                              readRDS("output/found.RatNr21_1000x_Th-20.rds"))

total.1000x <- cbind(unlist(total.1000x))[[1]]


venn.area1 <- length(unique(mcols(total.100x)$BC))
venn.area2 <- length(unique(mcols(total.1000x)$BC))

venn.n12 <- length(intersect(unique(mcols(total.100x)$BC),unique(mcols(total.1000x)$BC)))

venn.colors <- c("cornflower blue", "red") 
grid.newpage()
venn.plot <- draw.pairwise.venn(area1    = venn.area1,
                              area2    = venn.area2,
                              cross.area      = venn.n12,
                              scaled   = TRUE,
                              fill     = venn.colors,
                              alpha    = 0.3,
                              lty      = "blank",
                              cex      = 2,
                              cat.cex  = 2,
                              cat.col  = venn.colors,
                              category = c("100x", "1000x"))
grid.draw(venn.plot)


venn.area1 <- length(unique(names(total.100x)))
venn.area2 <- length(unique(names(total.1000x)))

venn.n12 <- length(intersect(unique(names(total.100x)),unique(names(total.1000x))))


grid.newpage()
venn.plot <- draw.pairwise.venn(area1    = venn.area1,
                                area2    = venn.area2,
                                cross.area      = venn.n12,
                                scaled   = TRUE,
                                fill     = venn.colors,
                                alpha    = 0.3,
                                lty      = "blank",
                                cex      = 2,
                                cat.cex  = 2,
                                cat.col  = venn.colors,
                                category = c("100x", "1000x"))
grid.draw(venn.plot)

lib.BC.counts <- data.table(mcols(complete.ranges)$tCount,mcols(complete.ranges)$BC)
lib.BC.counts <- lib.BC.counts[match(unique(lib.BC.counts$V2),lib.BC.counts$V2),]
lib.BC.counts
setorder(lib.BC.counts, V1)
lib.BC.counts$V3 <- 1:nrow(lib.BC.counts)
setkey(lib.BC.counts, V1)
lib.BC.counts.low <- lib.BC.counts[,ReadCount=min(V1), Position=min(V3), by=list(V1)] 
countPlot <- ggplot(lib.BC.counts,aes(x=V3, y=V1)) + geom_area() + 

  geom_histogram(bin=1, stat="identity")
, aes(fill = Library,y=ReadCount)

+theme_bw()+
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  facet_grid(GeneName~., scales = "free_x", space = "free_x") 
