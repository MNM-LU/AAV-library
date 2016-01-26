
#' ---
#' title: "Top 10 heatmap analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This is the final script presenting top 10 candidates as heatmap plots.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
opts_chunk$set(fig.width = 8, fig.height = 10.2)
opts_chunk$set(comment = NA)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(devtools))


#'Selection of relevant samples
#'===================
select.samples <- readRDS("data/allSamplesDataTable.RDS")
select.samples$Group[select.samples$Group== "293T_1000x"] <- "H293T_1000x"
select.samples$Group[select.samples$Group== "293T_100x"] <- "H293T_100x"
select.samples <- select.samples[-grep("4wks|PrimN_1000x_RNA",select.samples$Group),]

select.samples.merge <- data.table::copy(select.samples)
select.samples.merge$Lib <- "100x"
select.samples.merge$Lib[grep("1000x",select.samples.merge$Group)] <- "1000x"
select.samples.merge[,Group:=gsub("100x|1000x","",Group)]


#Connect highest scoring fragments
setkey(select.samples.merge,Group)
select.samples.merge.binPos <- select.samples.merge[c("CNS_Str","CNS_Th","CNS_Ctx","CNS_SN")]


setorder(select.samples.merge.binPos,Group,GeneName,start,width)
winWidth=1
windowTable <- select.samples.merge.binPos[,c("GeneName","start","width"), with = FALSE]
windowTable <- unique(windowTable, by=c("GeneName","start","width"))
windowTable <- windowTable[,(seq(width-winWidth+1)+start-1),by=c("GeneName","start","width")]
setnames(windowTable,"V1","winStart")
windowTable[,winEnd:=winStart+winWidth-1]
setkeyv(windowTable,c("GeneName","start","width"))
setkeyv(select.samples.merge.binPos,c("GeneName","start","width"))
select.samples.windowBin <- select.samples.merge.binPos[windowTable, allow.cartesian=TRUE]

setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))

select.samples.windowBin <- select.samples.windowBin[, list(Overlaps=.N,
                                                            seqlength=min(seqlength),
                                                            BCcount=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                                            NormCount=mean(log2(RNAcount+1)),
                                                            AnimalCount=length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                                            LUTnrs=paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                                            mainStruct=paste(unique(structure), collapse=","),
                                                            mainLibs=paste(unique(Lib), collapse=","),
                                                            libCount=length(unique(Lib)),
                                                            mismatches=median(mismatches)
), by=c("Group","GeneName","winStart","winEnd")]

select.samples.windowBin[,Score:=BCcount+AnimalCount+libCount]
setorder(select.samples.windowBin, Group,GeneName,winStart,-Score)



windowTable <- select.samples.windowBin[,c("Group","GeneName","seqlength"), with = FALSE]
windowTable <- windowTable[,seq(seqlength),by=c("Group","GeneName")]
setnames(windowTable,"V1","winStart")
setkeyv(windowTable,c("Group","GeneName","winStart"))
setkeyv(select.samples.windowBin,c("Group","GeneName","winStart"))
select.samples.windowBin.full <- data.table::copy(select.samples.windowBin)
select.samples.windowBin <- select.samples.windowBin[windowTable]
select.samples.windowBin$Score[is.na(select.samples.windowBin$Score)] <- 0

select.samples.windowBin$LUTnrs[is.na(select.samples.windowBin$LUTnrs)] <- seq(length(which(is.na(select.samples.windowBin$LUTnrs))))

select.samples.windowBin.unique <- unique(select.samples.windowBin, by=c("Group","LUTnrs","Score"))

setorder(select.samples.windowBin.unique, Group,GeneName,winStart)


local_max <- function(x) {
  which(diff(sign(diff(x)))==-2)+1}

select.samples.windowBin.locMax <- select.samples.windowBin.unique[local_max(select.samples.windowBin.unique$Score),]
setorder(select.samples.windowBin.locMax, Group,-Score,GeneName,winStart)
select.samples.windowBin.locMax <- unique(select.samples.windowBin.locMax, by=c("Group","LUTnrs"))

#'Make bins
#'===================


#Connect highest scoring fragments

setorder(select.samples.windowBin.locMax,Group,GeneName,winStart)
winWidth=7
windowTable <- select.samples.windowBin.locMax[,c("GeneName","winStart"), with = FALSE]
windowTable <- unique(windowTable, by=c("GeneName","winStart"))
windowTable <- windowTable[,seq((winStart-winWidth),(winStart+winWidth)),by=c("GeneName","winStart")]
setnames(windowTable,"V1","binBaseStart")
windowTable[,binBaseEnd:=binBaseStart+(2*winWidth)-1]
scoreSelect <- select.samples.windowBin.locMax[,c("Group","GeneName","winStart","Score"), with = FALSE]
setkeyv(windowTable,c("GeneName","winStart"))
setkeyv(scoreSelect,c("GeneName","winStart"))
scoreSelectBin <- scoreSelect[windowTable,allow.cartesian=TRUE]

scoreSelectBin[,mCount:=.N,by=c("GeneName","binBaseStart")]
scoreSelectBin[,oCount:=.N,by=c("GeneName","winStart","binBaseStart")]
scoreSelectBin[,offset:=abs(binBaseStart+winWidth-winStart)]
#scoreSelectBin <- scoreSelectBin[mCount>oCount,]
setorder(scoreSelectBin,GeneName,winStart,-mCount,-offset)
scoreSelectBinTop <- scoreSelectBin[, head(.SD, 1), by=c("GeneName","winStart")]
scoreSelect <- scoreSelectBinTop[,c("GeneName","winStart","binBaseStart","binBaseEnd"), with = FALSE]
setkeyv(scoreSelect,c("GeneName","winStart"))
setkeyv(select.samples.windowBin.locMax,c("GeneName","winStart"))
select.samples.windowBin.locMax.bin <- select.samples.windowBin.locMax[scoreSelect, allow.cartesian=TRUE]
setorder(select.samples.windowBin.locMax.bin,GeneName,-binBaseStart,-Score,Group)
select.samples.windowBin.locMerge <- select.samples.windowBin.locMax.bin[, head(.SD, 1), by=c("GeneName","binBaseStart","Group")]


#'Selection of top twenty fragments per sample
#'===================


setorder(select.samples.windowBin.locMax,Group,-Score,-AnimalCount,-BCcount,-libCount,-NormCount,GeneName,winStart)
setkey(select.samples.windowBin.locMax,Group)

select.samples.topTwenty <- select.samples.windowBin.locMax[, head(.SD, 10), by=Group]
#'Selectonly the active samples
# setkey(select.samples.topTwenty,Group)
# select.samples.topTwenty <- select.samples.topTwenty[c("CNS_Th","CNS_Ctx","CNS_SN")]

select.samples.topTwenty <- select.samples.topTwenty[,c("GeneName","winStart"),with=FALSE]
select.samples.topTwenty <- unique(select.samples.topTwenty, by=c("GeneName","winStart"))


setkeyv(select.samples.topTwenty,c("GeneName","winStart"))
setkeyv(select.samples.windowBin.full,c("GeneName","winStart"))
select.samples.windowBin.allTop <- select.samples.windowBin.full[select.samples.topTwenty, nomatch=0]
select.samples.windowBin.allTop <- unique(select.samples.windowBin.allTop, by=c("Group","LUTnrs","Score"))

select.samples.windowBin.allTop[,GeneAA:=paste(GeneName," [",winStart,"]", sep="")]
setorder(select.samples.windowBin.allTop,Group,-Score)

#'Plotting ranked order
#'===================

make_bipartite_graph(select.samples.windowBin.allTop$Group, edges, directed = FALSE)
setkey(select.samples.windowBin.allTop,Group)

v1 <- select.samples.windowBin.allTop["CNS_Str"]$GeneAA
v2 <- select.samples.windowBin.allTop["CNS_Th"]$GeneAA
v3 <- select.samples.windowBin.allTop["CNS_Ctx"]$GeneAA
v4 <- select.samples.windowBin.allTop["CNS_SN"]$GeneAA


o <- 0.05
DF <- data.table(x = c(rep(1, length(v1)), rep(2, length(v2)),rep(3, length(v3)), rep(4, length(v4))),
                 x1 = c(rep(1, length(v1)),rep(2, length(v2)), rep(3, length(v3)), rep(4, length(v4))),
                 y = c(rev(seq_along(v1)), rev(seq_along(v2)), rev(seq_along(v3)), rev(seq_along(v4))),
                 g = c(v1, v2, v3, v4))
DF[,groupMax:=max(y), by=x]
allMax <- max(DF$y)
DF[,y:=y-groupMax+allMax]


library(ggplot2)
library(grid)
ggplot(DF, aes(x=x, y=y, group=g, label=g)) +
  geom_path(aes(x=x1), 
            size=0.5, color="blue") +
  geom_text(size=2.6) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

