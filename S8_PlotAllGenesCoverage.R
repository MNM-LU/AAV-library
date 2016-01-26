
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
#+ setup, include=FALSE

opts_chunk$set(fig.width = 8, fig.height = 10.2) #Full height 11
opts_chunk$set(comment = NA)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))

#'Generation of infective library
#'===================
all.samples <- readRDS("data/allSamplesDataTable.RDS")

all.samples$Group[all.samples$Group== "293T_1000x"] <- "H293T_1000x"
all.samples$Group[all.samples$Group== "293T_100x"] <- "H293T_100x"

#'Plotting function
#'===================

plotPair <- function(topSample,bottomSample,size.bin=1,winWidth=1,NormalizePlot=TRUE, PlotBC=TRUE) {
# Select samples
#===================

#   topSample <- "CNS1000x_Ctx"
#   bottomSample <- "CNS100x_Ctx"
#   filterBC <- FALSE
#   filterAnimal <- FALSE
#   AnimaladjustPlot <- FALSE
#   NormalizePlot <- TRUE
#   size.bin <- 1
#   winWidth=1
#   PlotBC=TRUE
  
fill.values <- eval(parse(text=paste("c(", topSample,"= rgb(38,64,135, maxColorValue = 255), ",
                                     bottomSample,"= rgb(157,190,217, maxColorValue = 255))",sep="")))
setkey(all.samples,Group)
select.samples <- all.samples[J(names(fill.values))] #Select the two compared groups
select.samples[,RNAcount:=log2(RNAcount+1)]

# if (PlotBC){
# select.samples[,c("meanCount","SDcount","minCount"):=list(mean(RNAcount),
#                                                sd(RNAcount),
#                                                min(RNAcount)),by="Group"]
# select.samples <- select.samples[RNAcount>=minCount+(2*SDcount),] #-SDcount
# select.samples[,c("meanCount","SDcount","minCount"):=NULL]
# }




if (winWidth > 0) {
  setorder(select.samples,Group,GeneName,start,width)
  
  windowTable <- select.samples[,c("GeneName","start","width"), with = FALSE]
  windowTable <- unique(windowTable, by=c("GeneName","start","width"))
  windowTable <- windowTable[,(seq(width-winWidth+1)+start-1),by=c("GeneName","start","width")]
  setnames(windowTable,"V1","winStart")
  windowTable[,winEnd:=winStart+winWidth-1]
  setkeyv(windowTable,c("GeneName","start","width"))
  setkeyv(select.samples,c("GeneName","start","width"))
  select.samples.windowBin <- select.samples[windowTable, allow.cartesian=TRUE]
  select.samples.windowBin[,AAproc:=winStart/seqlength*100]
  
setkey(select.samples.windowBin,Group)
select.samples.windowBin <- select.samples.windowBin[J(names(fill.values))] #Select the two compared groups
setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))
select.samples.windowBin <- select.samples.windowBin[, list(Overlaps=.N,
                                        seqlength=min(seqlength),
                                        AAproc = min(AAproc),
                                        BC = paste(t(BC), collapse=","),
                                        Animals = paste(t(Animals), collapse=","),
                                        LUTnrs = paste(t(LUTnrs), collapse=","),
                                        RNAcount = sum(RNAcount)
), by=c("Group","GeneName","winStart","winEnd")]

plot.data.dt <- unique(select.samples.windowBin, by=c("Group","GeneName","winStart","winEnd"))

} else {
  plot.data.dt <- data.table::copy(select.samples)
}

#===================
#Binning of data
#===================
FullLength <- 100
position <- seq(0,FullLength,size.bin)
plot.data.dt[,bin:=findInterval(AAproc, position)]

plot.data.bin <- plot.data.dt[, list(.N,seqlength=min(seqlength),
                                     AAproc = position[findInterval(mean(AAproc),position)],
                                     BCsum=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                     AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                     LUTnrs = paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                     NormCount = sum(RNAcount)/seqlength*FullLength
), by=c("Group","GeneName","bin")]
plot.data.bin <- unique(plot.data.bin, by=c("Group","GeneName","bin"))

plot.data.bin[,BCanim:=as.double(BCsum*AnimalCount)]



#===================
#Filtration parameters
#===================

if (NormalizePlot) {
  plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount)
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount)
}


if (PlotBC && NormalizePlot) {
  plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$BCanim)
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim)
}


plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount*-1 #This line flips the values for the second group

plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$BCanim*-1 #This line flips the values for the second group

#===================
#Output plot
#===================

if (PlotBC) {outVar <- "BCanim"} else {outVar <- "NormCount"}

plot.out <- eval(parse(text=paste("ggplot(plot.data.bin,aes(x=AAproc,y=",outVar,", fill = Group))",sep="")))
plot.out <- plot.out + 
  geom_bar(stat="identity", position="identity")+theme_bw() +
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  scale_x_continuous(limit=c(0,100), breaks=c(seq(0,100,20)), expand =c(0,0)) +
  facet_wrap(~ GeneName, ncol=5)+   
  theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),
        legend.position="bottom",
        legend.margin=unit(-0.5,"cm"),
        legend.key.height=unit(0, "cm"),
        plot.background=element_rect(fill="white"),
        axis.text = element_text(size = rel(0.45)),
        axis.ticks = element_line(size = rel(0.5)),
        axis.ticks.length = unit(.05, "cm"),
        strip.text.x = element_text(size = rel(0.5), colour = "black", 
                                    angle = 0, lineheight=0.1, vjust=0.1),
        strip.background = element_blank(),
        strip.background = element_rect(size = 0),
        panel.margin.y = unit(-0.15, "cm"),
        panel.margin.x = unit(0, "cm"))


#===================
# Sort and select top samples
#===================

select.samples.binPos <<- select.samples
setkeyv(select.samples.binPos,c("Group","structure","Sequence"))
setorder(select.samples.binPos,Group,structure,Sequence,GeneName)
select.samples.binPos <- unique(select.samples.binPos, by=c("Group","structure","Sequence")) 
#Due to key, this removes replicates if identical sequence mapped to multiple genes

setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))
select.samples.binPos[,c("BCcount","NormCount","AnimalCount","LUTnrs","mainStruct","mismatches"):=
                        list(length(table(strsplit(paste(t(BC), collapse=","), ","))),
                             sum(NormCount),
                             length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                             paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                             paste(unique(structure), collapse=","),
                             median(mismatches)), by=key(select.samples.binPos)]

select.samples.binPos <- unique(select.samples.binPos, by=c("Group","NormCount","LUTnrs"))
select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount",
                                                  "BCcount","AnimalCount","LUTnrs","mainStruct",
                                                  "mismatches"), with = FALSE]

if (PlotBC){
  select.samples.binPos[,BCanim:=BCcount*AnimalCount]
  setorder(select.samples.binPos,Group,-BCanim,-BCcount,-AnimalCount,-NormCount)
} else {
setorder(select.samples.binPos,Group,-NormCount,-BCcount,-AnimalCount)
}
setkey(select.samples.binPos,Group)
select.samples.top <- select.samples.binPos[, head(.SD, 20), by=Group]
top.sample <- select.samples.top[J(names(fill.values)[1])]
bottom.sample <- select.samples.top[J(names(fill.values)[2])]
top.sample[,c("Group"):=NULL]
setnames(top.sample, "GeneName", names(fill.values)[1])
bottom.sample[,c("Group"):=NULL]
setnames(bottom.sample, "GeneName", names(fill.values)[2])

out.list <- list(plot=plot.out,
                 plotBin=plot.data.bin,
                 top=top.sample,
                 bottom=bottom.sample)

return(out.list)
}


#'Analyze samples
#'===================

plotPair("infectiveLib","totalLib")$plot
#' 100x analysis
out.plot.list <- plotPair("CNS100x_Trsp","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS100x_Th","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS100x_Ctx","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS100x_SN","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS100x_SN","CNS100x_Th")$plot
plotPair("CNS100x_Ctx","CNS100x_Th")$plot
plotPair("CNS100x_SN","CNS100x_Ctx")$plot

#' 1000x analysis
out.plot.list <- plotPair("CNS1000x_Trsp","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS1000x_Th","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS1000x_Ctx","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS1000x_SN","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS1000x_SN","CNS1000x_Th")$plot
plotPair("CNS1000x_Ctx","CNS1000x_Th")$plot
plotPair("CNS1000x_SN","CNS1000x_Ctx")$plot

#' 1000x vs 100x analysis
plotPair("CNS1000x_Trsp","CNS100x_Trsp")$plot
plotPair("CNS1000x_Str","CNS100x_Str")$plot
plotPair("CNS1000x_Th","CNS100x_Th")$plot
plotPair("CNS1000x_Ctx","CNS100x_Ctx")$plot
plotPair("CNS1000x_SN","CNS100x_SN")$plot

#' 8wks vs 4wks analysis

plotPair("CNS100x_Str_4wks","CNS100x_Str")$plot
plotPair("CNS100x_Th_4wks","CNS100x_Th")$plot
plotPair("CNS100x_Ctx_4wks","CNS100x_Ctx")$plot
plotPair("CNS100x_SN_4wks","CNS100x_SN")$plot
plotPair("PerN100x_SC_4wks","PerN100x_Mu_4wks")$plot

plotPair("CNS1000x_Str_4wks","CNS1000x_Str")$plot
plotPair("CNS1000x_Th_4wks","CNS1000x_Th")$plot
plotPair("CNS1000x_Ctx_4wks","CNS1000x_Ctx")$plot
plotPair("CNS1000x_SN_4wks","CNS1000x_SN")$plot

#' In vitro analysis

out.plot.list <- plotPair("PrimN_1000x","PrimN_100x")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("H293T_1000x","H293T_100x")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")



#' Old binning analysis

plotPair("infectiveLib","totalLib",PlotBC=FALSE)$plot
#' 100x analysis
out.plot.list <- plotPair("CNS100x_Trsp","CNS100x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS100x_Th","CNS100x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS100x_Ctx","CNS100x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS100x_SN","CNS100x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS100x_SN","CNS100x_Th",PlotBC=FALSE)$plot
plotPair("CNS100x_Ctx","CNS100x_Th",PlotBC=FALSE)$plot
plotPair("CNS100x_SN","CNS100x_Ctx",PlotBC=FALSE)$plot

#' 1000x analysis
out.plot.list <- plotPair("CNS1000x_Trsp","CNS1000x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS1000x_Th","CNS1000x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS1000x_Ctx","CNS1000x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS1000x_SN","CNS1000x_Str",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS1000x_SN","CNS1000x_Th",PlotBC=FALSE)$plot
plotPair("CNS1000x_Ctx","CNS1000x_Th",PlotBC=FALSE)$plot
plotPair("CNS1000x_SN","CNS1000x_Ctx",PlotBC=FALSE)$plot

#' 1000x vs 100x analysis
plotPair("CNS1000x_Trsp","CNS100x_Trsp",PlotBC=FALSE)$plot
plotPair("CNS1000x_Str","CNS100x_Str",PlotBC=FALSE)$plot
plotPair("CNS1000x_Th","CNS100x_Th",PlotBC=FALSE)$plot
plotPair("CNS1000x_Ctx","CNS100x_Ctx",PlotBC=FALSE)$plot
plotPair("CNS1000x_SN","CNS100x_SN",PlotBC=FALSE)$plot

#' 8wks vs 4wks analysis

plotPair("CNS100x_Str_4wks","CNS100x_Str",PlotBC=FALSE)$plot
plotPair("CNS100x_Th_4wks","CNS100x_Th",PlotBC=FALSE)$plot
plotPair("CNS100x_Ctx_4wks","CNS100x_Ctx",PlotBC=FALSE)$plot
plotPair("CNS100x_SN_4wks","CNS100x_SN",PlotBC=FALSE)$plot
plotPair("PerN100x_SC_4wks","PerN100x_Mu_4wks",PlotBC=FALSE)$plot

plotPair("CNS1000x_Str_4wks","CNS1000x_Str",PlotBC=FALSE)$plot
plotPair("CNS1000x_Th_4wks","CNS1000x_Th",PlotBC=FALSE)$plot
plotPair("CNS1000x_Ctx_4wks","CNS1000x_Ctx",PlotBC=FALSE)$plot
plotPair("CNS1000x_SN_4wks","CNS1000x_SN",PlotBC=FALSE)$plot

#' In vitro analysis

out.plot.list <- plotPair("PrimN_1000x","PrimN_100x",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("H293T_1000x","H293T_100x",PlotBC=FALSE)
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

plotPair("infectiveLib","totalLib",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Str","CNS100x_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Th","CNS100x_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Ctx","CNS100x_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_SN","CNS100x_SN",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("CNS100x_SN","CNS100x_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS100x_Ctx","CNS100x_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS100x_Th","CNS100x_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("CNS1000x_SN","CNS1000x_Ctx",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Ctx","CNS1000x_Th",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Th","CNS1000x_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("CNS100x_Trsp","CNS100x_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Trsp","CNS1000x_Str",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("CNS1000x_Trsp","CNS100x_Trsp",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

plotPair("PrimN_1000x","PrimN_100x",size.bin=2,winWidth=0,PlotBC=FALSE)$plot
plotPair("H293T_1000x","H293T_100x",size.bin=2,winWidth=0,PlotBC=FALSE)$plot

devtools::session_info()
