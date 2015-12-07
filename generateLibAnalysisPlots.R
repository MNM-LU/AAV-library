suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(data.table))

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


fill.values <- c("Array" = rgb(193,210,234, maxColorValue = 255), "Lib" = rgb(157,190,217, maxColorValue = 255), "AAV" = rgb(38,64,135, maxColorValue = 255), "Str" = rgb(157,190,217, maxColorValue = 255),"SN" = rgb(193,210,234, maxColorValue = 255), "None" = rgb(255,255,255, maxColorValue = 255))


ggplot(output.table) + 
  geom_rect(aes(fill=NameSN, ymax=SNEnd, ymin=SNStart, xmax=13, xmin=11)) +
  geom_rect(aes(fill=NameStr, ymax=StrEnd, ymin=StrStart, xmax=11, xmin=9)) +
  geom_rect(aes(fill=NameAAV, ymax=AAVend, ymin=AAVStart, xmax=9, xmin=7)) +
  geom_rect(aes(fill=NameLib, ymax=LibEnd, ymin=LibStart, xmax=7, xmin=4)) +
  geom_rect(aes(fill=NameArray, ymax=ArrayEnd, ymin=ArrayStart, xmax=4, xmin=0)) +
  xlim(c(0, 13)) + 
  scale_fill_manual(name = "Library", values = fill.values) +
  theme(aspect.ratio=1) 





grid.newpage()
venn.plot <- draw.triple.venn(area1    = venn.area1,
                              area2    = venn.area2,
                              area3    = venn.area3,
                              n12      = venn.n12,
                              n23      = venn.n23,
                              n13      = venn.n13,
                              n123     = venn.n123,
                              scaled   = FALSE,
                              fill     = c("blue", "red", "green"),
                              alpha    = 0.3,
                              lty      = "blank", ,
                              cex      = 2,
                              cat.cex  = 2,
                              cat.col  = c("blue", "red", "green"),
                              category = c("CustomArry", "PlasmidLib", "InfectiveAAV"),
                              cat.pos  = c(0, 40, 250),
                              cat.dist = c(0.05, 0.05, 0.05))
grid.draw(venn.plot)


ggplot(browsers) + 
  geom_rect(aes(fill=version, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(aes(fill=browser, ymax=ymax, ymin=ymin, xmax=3, xmin=0)) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  coord_polar(theta="y")  


