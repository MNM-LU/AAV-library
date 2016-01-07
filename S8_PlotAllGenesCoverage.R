
#' ---
#' title: "Two sample analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This is the final script presenting top candidates and overview plots.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
opts_chunk$set(fig.width = 8, fig.height = 11) #Full height 11
opts_chunk$set(comment = NA)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(devtools))
library(grid)
# Colors
# "springgreen4"
#rgb(38,64,135, maxColorValue = 255)
#rgb(157,190,217, maxColorValue = 255)

#'Setup parameters
#'===================
all.samples <- readRDS("data/normalizedSampleRanges.RDS")
fill.values <- c("CNS100x_SN" = rgb(38,64,135, maxColorValue = 255), 
                 "CNS100x_Str" = rgb(157,190,217, maxColorValue = 255))
filterBC <- FALSE
filterAnimal <- FALSE
AnimaladjustPlot <- FALSE
NormalizePlot <- TRUE


#'Generation of infective library
#'===================
total.AAV.samples <- all.samples[!(mcols(all.samples)$Group %in% "totalLib")]
total.AAV.samples <- total.AAV.samples[-grep("4wks",mcols(total.AAV.samples)$Group)]
mcols(total.AAV.samples)$Group <- "infectiveLib"
all.samples <- append(all.samples,total.AAV.samples)
select.samples <- all.samples[mcols(all.samples)$Group %in% names(fill.values)]
trim.names <- data.table(GeneName=levels(seqnames(select.samples)))
trim.names[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
trim.names$GeneName <- gsub("/","-",trim.names$GeneName)
trim.names$GeneName <- gsub("_","-",trim.names$GeneName)
levels(seqnames(select.samples)) <- trim.names$GeneName
select.samples <- select.samples[mcols(select.samples)$LV < 3 & mcols(select.samples)$mCount/mcols(select.samples)$tCount > 0.5]

geneTable <- data.frame(seqlengths(select.samples))
colnames(geneTable) <- "SeqLength"
geneTable$GeneName <- row.names(geneTable)

plot.data <- data.frame(GeneName=seqnames(select.samples),
                        Group=factor(mcols(select.samples)$Group, levels = names(fill.values)),
                        AApos=start(select.samples)+(width(select.samples)/2),
                        Animals=mcols(select.samples)$Animals,
                        BC=mcols(select.samples)$BC,
                        RNAcount=mcols(select.samples)$RNAcount,
                        NormCount=mcols(select.samples)$NormCount,
                        stringsAsFactors = FALSE)
plot.data <- merge(plot.data,geneTable)
plot.data$AAproc <- (plot.data$AApos/plot.data$SeqLength)*100

#'Binning of data
#'===================

size.bin <- 2
FullLength <- 100
position <- seq(0,FullLength,size.bin)
plot.data.dt <- data.table(plot.data)
plot.data.dt[,bin:=findInterval(AAproc, position)]
plot.data.bin <- plot.data.dt[, list(.N,SeqLength=min(SeqLength),
                                        AAproc = position[findInterval(min(AAproc),position)],
                                        BCmean=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                        AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                        NormCount = log2((sum(RNAcount)/SeqLength*FullLength)+1)
                                        ), by=c("bin","GeneName","Group")]



knitr::kable(plot.data.bin[1:10,], format = "markdown")
if (NormalizePlot) {
  plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount)
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount)
}

if (filterBC) {
  plot.data.bin <- plot.data.bin[plot.data.bin$BCmean > 1,]
  select.samples <- select.samples[mcols(select.samples)$BC>1]
}

if (filterAnimal) {
  plot.data.bin <- plot.data.bin[plot.data.bin$AnimalCount > 1,]
}

if (AnimaladjustPlot) {
  plot.data.bin$NormCount <- plot.data.bin$NormCount*plot.data.bin$AnimalCount/plot.data.bin$SeqLength*FullLength
}

plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount*-1


ggplot(plot.data.bin,aes(x=AAproc,y=NormCount, fill = Group))+geom_bar(stat="identity", position="identity")+theme_bw()+
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  scale_x_continuous(limit=c(0,100), breaks=c(seq(0,100,20)), expand =c(0,0)) +
  facet_wrap(~ GeneName, ncol=5)+   
  theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),
    legend.position="bottom",
    legend.margin=unit(-0.5,"cm"),
    legend.key.height=unit(0, "cm"),
    plot.background=element_rect(fill="white"),
    axis.text = element_text(size = rel(0.5)),
    axis.ticks = element_line(size = rel(0.5)),
    axis.ticks.length = unit(.05, "cm"),
    strip.text.x = element_text(size = rel(0.5), colour = "black", angle = 0, lineheight=0.1, vjust=0.1),
    strip.background = element_blank(),
    strip.background = element_rect(size = 0),
    panel.margin.y = unit(-0.15, "cm"),
    panel.margin.x = unit(0, "cm"))


  #facet_grid(GeneName~., scales = "free_x", space = "free_x") 

select.samples.gr <- granges(select.samples)
mcols(select.samples.gr) <- cbind(mcols(select.samples)[,c(1,2,3,4,5,6,9)],
                                  DataFrame(Animals=unlist(lapply(strsplit(mcols(select.samples)$Animals, ","),function(x)length(table(x))))))

o = order(-mcols(select.samples.gr)$NormCount)
select.samples.gr <- select.samples.gr[o]
top.sample <- select.samples.gr[mcols(select.samples.gr)$Group %in% names(fill.values)[1]]
bottom.sample <- select.samples.gr[mcols(select.samples.gr)$Group %in% names(fill.values)[2]]
#+ setup, fontsize: 7pt
if (length(top.sample) >=1){
out <- top.sample
out.2 <- data.table(Fragment=names(out),
                    SeqCount=0,
                  Gene=as.character(seqnames(out)),
                  startAA=(start(out)+2)/3)
out.2 <- data.frame(out.2,as.data.frame(mcols(out)))
out.2$BC <- unlist(lapply(strsplit(out.2$BC, ","),function(x) length(unique(x))))
out.2[match(names(table(out.2$Fragment)),out.2$Fragment),"SeqCount"] <- as.integer(table(out.2$Fragment))
out.2 <- out.2[out.2$SeqCount>=1,]
return(knitr::kable(out.2[1:min(15,nrow(out.2)),c(1,3,6)], format = "markdown"))
}
if (length(top.sample) >=1){
  return(knitr::kable(out.2[1:min(15,nrow(out.2)),c(-1,-6)], format = "markdown"))
}

if (length(bottom.sample) >=1) {
  out <- bottom.sample
  out.2 <- data.table(Fragment=names(out),
                      SeqCount=0,
                      Gene=as.character(seqnames(out)),
                      startAA=(start(out)+2)/3)
  out.2 <- data.frame(out.2,as.data.frame(mcols(out)))
  out.2$BC <- unlist(lapply(strsplit(out.2$BC, ","),function(x) length(unique(x))))
  out.2[match(names(table(out.2$Fragment)),out.2$Fragment),"SeqCount"] <- as.integer(table(out.2$Fragment))
  out.2 <- out.2[out.2$SeqCount>=1,]
  return(knitr::kable(out.2[1:min(15,nrow(out.2)),c(1,3,6)], format = "markdown"))
}
if (length(bottom.sample) >=1) {
  return(knitr::kable(out.2[1:min(15,nrow(out.2)),c(-1,-6)], format = "markdown"))
}

devtools::session_info()