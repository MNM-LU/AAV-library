
#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This is the final script presenting top candidates and overview plots.  
suppressPackageStartupMessages(library(knitr))

opts_chunk$set(fig.width = 5, fig.height = 5) #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)
#+ setup, include=FALSE
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(beanplot))
suppressPackageStartupMessages(library(scales))

#'Generation plots for library purity
#'===================
complete.ranges <- readRDS("output/completeLibraryRanges.rds")
purity.table <- data.table(mcols(complete.ranges)$mCount/mcols(complete.ranges)$tCount)
purity.table$BCwidth <- width(mcols(complete.ranges)$BC)

beanplot(data = purity.table$V1, ll = 0.04, what=c(1,1,1,0), bw = "nrd0", log="", 
         main = "Best end recombination analysis", ylab = "Recombination fraction [1 = no recombination]", 
         border = NA, col = list("black", c("grey", "white")))
legend("bottomleft", fill = c("black"), legend = c("Best end"))

fill.values <- c("A" = rgb(193,210,234, maxColorValue = 255), "B" = rgb(157,190,217, maxColorValue = 255), "c" = rgb(38,64,135, maxColorValue = 255), "D" = rgb(157,190,217, maxColorValue = 255),"E" = rgb(193,210,234, maxColorValue = 255))
names(fill.values) <- c("18","19","20","21","22")

ggplot(purity.table,aes(x = BCwidth, fill = as.character(BCwidth)))+geom_histogram(binwidth=1, boundary = 0.5)+theme_bw() +
  scale_fill_manual(name = "Width", values = fill.values) +
  scale_x_continuous(limit=c(17,23), breaks=c(seq(18,22,1)), expand =c(0,0))

opts_chunk$set(fig.width = 8, fig.height = 8)
#'Plot Venn diagrams of fragments
#'===================
load("data/LUTdna.rda")
complete.library <- readRDS("data/allSamplesDataTable.RDS")
setkey(complete.library,Group)
#complete.library <- complete.library[-grep("4wk",Group)]
seq.arry <- LUT.dna$LUTnr
seq.lib <- unique(complete.library[J("totalLib")]$LUTnr)
seq.AAV <- unique(complete.library[J("infectiveLib")]$LUTnr)
seq.str <- unique(complete.library[grep("Str",Group)]$LUTnr)
seq.Trsp <- unique(complete.library[grep("SN|Ctx|Th",Group)]$LUTnr)

venn.area1 <- length(seq.arry)
venn.area2 <- length(seq.lib)
venn.area3 <- length(seq.AAV)
venn.area3 <- length(seq.str)
venn.area3 <- length(seq.Trsp)


isect.Str_Trsp <- length(intersect(seq.str, seq.Trsp))

venn.n12 <- length(intersect(seq.arry,seq.lib))
venn.n23 <- length(intersect(seq.lib,seq.AAV))
venn.n13 <- length(intersect(seq.arry,seq.AAV))
venn.n123 <- length(intersect(intersect(seq.arry,seq.lib),seq.AAV))


output.table <- data.frame(NameArray=character(),
                           NameLib=character(),
                           NameAAV=character(),
                           NameStr=character(),
                           NameTrsp=character(),
                           ArrayStart=numeric(), 
                           ArrayEnd=numeric(), 
                           LibStart=numeric(), 
                           LibEnd=numeric(), 
                           AAVStart=numeric(),
                           AAVend=numeric(),
                           StrStart=numeric(),
                           StrEnd=numeric(),
                           TrspStart=numeric(),
                           TrspEnd=numeric(),
                           stringsAsFactors=FALSE) 
output.table[1:3,1] <- c("Array","None","None")
output.table[1:3,2] <- c("Lib","None","None")
output.table[1:3,3] <- c("AAV","None","None")
output.table[1:3,4] <- c("Str","None","None")
output.table[1:3,5] <- c("None","Trsp","None")
output.table[1:3,6:15] <- 0
output.table$ArrayEnd[1] <- output.table$LibEnd[2] <- output.table$AAVend[2] <- output.table$StrEnd[2] <- output.table$TrspEnd[3] <- length(seq.arry)
output.table$LibStart[2] <- output.table$LibEnd[1] <- length(intersect(seq.arry,seq.lib))
output.table$AAVStart[2] <- output.table$AAVend[1] <- length(intersect(seq.lib,seq.AAV))
output.table$StrStart[2] <- output.table$StrEnd[1] <- length(intersect(seq.AAV, seq.str))
output.table$TrspStart[2] <- output.table$TrspEnd[1] <- length(seq.str)-length(intersect(seq.str, seq.Trsp))
output.table$TrspStart[3] <- output.table$TrspEnd[2] <- output.table$TrspStart[2] + length(seq.Trsp)


fill.values <- c("Array" = rgb(193,210,234, maxColorValue = 255), "Lib" = rgb(157,190,217, maxColorValue = 255), "AAV" = rgb(38,64,135, maxColorValue = 255), "Str" = rgb(157,190,217, maxColorValue = 255),"Trsp" = rgb(193,210,234, maxColorValue = 255), "None" = rgb(255,255,255, maxColorValue = 255, alpha = 0))


ggplot(output.table) + 
  scale_x_continuous(limit=c(0,10), breaks=c(seq(1,9,2)), expand =c(0,0)) + 
  scale_fill_manual(name = "Library", values = fill.values) +
  theme(aspect.ratio=1) + 
  geom_rect(data=output.table, aes(fill=NameTrsp, ymax=TrspEnd, ymin=TrspStart, xmax=10, xmin=8)) +
  geom_rect(data=output.table, aes(fill=NameStr, ymax=StrEnd, ymin=StrStart, xmax=8, xmin=6)) +
  geom_rect(data=output.table, aes(fill=NameAAV, ymax=AAVend, ymin=AAVStart, xmax=6, xmin=4)) +
  geom_rect(data=output.table, aes(fill=NameLib, ymax=LibEnd, ymin=LibStart, xmax=4, xmin=2)) +
  geom_rect(data=output.table, aes(fill=NameArray, ymax=ArrayEnd, ymin=ArrayStart, xmax=2, xmin=0)) +
  coord_polar(theta="y") 

opts_chunk$set(fig.width = 5, fig.height = 5)
#'Barcode Venn diagrams for 100x and 1000x libraries
#'===================
total.100x <- unique(complete.library[grep("100x",Group)]$BC)
total.1000x <- unique(complete.library[grep("1000x",Group)]$BC)


venn.area1 <- length(total.100x)
venn.area2 <- length(total.1000x)

venn.n12 <- length(intersect(total.100x,total.1000x))

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

#'Fragment Venn diagrams for 100x and 1000x libraries
#'===================
total.100x <- unique(complete.library[grep("100x",Group)]$Sequence)
total.1000x <- unique(complete.library[grep("1000x",Group)]$Sequence)


venn.area1 <- length(total.100x)
venn.area2 <- length(total.1000x)

venn.n12 <- length(intersect(total.100x,total.1000x))


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
