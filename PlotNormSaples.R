suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))

all.samples <- readRDS("output/normalizedSampleRanges.RDS")


unique(mcols(all.samples)$Sample)


plot.samples <- all.samples[mcols(all.samples)$Sample == "CNS1000x_SN" |
                              mcols(all.samples)$Sample == "CNS1000x_Ctx" |
                              mcols(all.samples)$Sample == "PerN100x_SC",]


plotData <- data.table(data.frame(Sample= mcols(plot.samples)$Sample,
                       GeneName=seqnames(plot.samples),
                       SamplePos = start(plot.samples)+(width(plot.samples)/2),
                       Score=mcols(plot.samples)$NormCount,
                       stringsAsFactors = FALSE))
plotData[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]


fill.values <- c("A" = rgb(157,190,217, maxColorValue = 255), "B" = rgb(38,64,135, maxColorValue = 255), "C" = "springgreen4")
names(fill.values) <- unique(plotData$Sample)

ggplot(plotData,aes(x=SamplePos,y=Score))+geom_histogram(bin=1, stat="identity", aes(fill = Sample,y=Score))+theme_bw()+
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  facet_grid(GeneName~., scales = "free_x", space = "free_x") 

one.samples <- all.samples[mcols(all.samples)$Sample == "CNS1000x_Ctx",]
mcols(one.samples)$NormCount/log2(mcols(one.samples)$RNAcount+1)
o = order(-mcols(one.samples)$NormCount)
one.samples <- one.samples[o]
one.samples[1:10]