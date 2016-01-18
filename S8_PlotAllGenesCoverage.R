
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

opts_chunk$set(fig.width = 8, fig.height = 11) #Full height 11
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
all.samples <- readRDS("data/normalizedSampleRangesDefined.RDS")
total.AAV.samples <- all.samples[!(mcols(all.samples)$Group %in% "totalLib")]
total.AAV.samples <- total.AAV.samples[-grep("4wks",mcols(total.AAV.samples)$Group)]
mcols(total.AAV.samples)$Group <- "infectiveLib"
all.samples <- append(all.samples,total.AAV.samples)

mcols(all.samples)$Sequence <- names(all.samples)
names(all.samples) <- make.names(names(all.samples), unique=TRUE)
length.Table <- data.table(seqnames=names(seqlengths(all.samples)),
                           seqlength=seqlengths(all.samples), key="seqnames")
all.samples <- data.table(as.data.frame(all.samples), key="seqnames")
all.samples[,c("strand","qwidth","cigar","njunc","end"):=NULL]
all.samples <- all.samples[length.Table] #A data.table merge to match seqlengths to their respective seqnames
all.samples[, c("Category", "Protein", "Origin", 
               "Extra", "Number","GeneName") := tstrsplit(seqnames, ",", fixed=TRUE)]
all.samples[, c("seqnames","Protein", "Origin", 
                "Extra", "Number") := NULL]
all.samples[, GeneName := gsub("/|_","-",GeneName)]
all.samples[,start:=(start+2)/3]
all.samples[,width:=(width)/3]
all.samples[,seqlength:=seqlength/3]
all.samples[,AA:=start+(width/2)]
all.samples[,AAproc:=AA/seqlength*100]

saveRDS(all.samples,file="data/allSamplesDataTable.RDS")

all.samples$Group[all.samples$Group== "293T_1000x"] <- "H293T_1000x"
all.samples$Group[all.samples$Group== "293T_100x"] <- "H293T_100x"

#'Plotting function
#'===================

plotPair <- function(topSample,bottomSample,size.bin=2,winWidth=1,NormalizePlot=TRUE) {
# Select samples
#===================

#   topSample <- "CNS1000x_SN"
#   bottomSample <- "CNS100x_SN"
#   filterBC <- FALSE
#   filterAnimal <- FALSE
#   AnimaladjustPlot <- FALSE
#   NormalizePlot <- TRUE
#   size.bin <- 1
#   winWidth=5
  
fill.values <- eval(parse(text=paste("c(", topSample,"= rgb(38,64,135, maxColorValue = 255), ",
                                     bottomSample,"= rgb(157,190,217, maxColorValue = 255))",sep="")))
setkey(all.samples,Group)
select.samples <- all.samples[J(names(fill.values))] #Select the two compared groups

if (winWidth > 1) {
  setorder(select.samples,Group,GeneName,start,width)
  
  windowTable <- select.samples[,c("GeneName","start","width"), with = FALSE]
  windowTable <- unique(windowTable, by=c("GeneName","start","width"))
  windowTable <- windowTable[,(seq(width-winWidth)+start-1),by=c("GeneName","start","width")]
  setnames(windowTable,"V1","winStart")
  windowTable[,winEnd:=winStart+winWidth]
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
                                        RNAcount = sum(RNAcount)
), by=c("Group","GeneName","winStart","winEnd")]

plot.data.dt <- unique(select.samples.windowBin)

} else {
  plot.data.dt <- copy(select.samples)
}

#===================
#Binning of data
#===================
FullLength <- 100
position <- seq(0,FullLength,size.bin)
plot.data.dt[,bin:=findInterval(AAproc, position)]
plot.data.bin <- plot.data.dt[, list(.N,seqlength=min(seqlength),
                                     AAproc = position[findInterval(min(AAproc),position)],
                                     BCmean=unlist(lapply(strsplit(paste(BC, collapse=","), ","),function(x) length(unique(x)))),
                                     AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                     NormCount = log2((sum(RNAcount)/seqlength*FullLength)+1)
), by=c("Group","GeneName","bin")]
plot.data.bin <- unique(plot.data.bin, by=c("Group","GeneName","bin"))

#===================
#Filtration parameters
#===================

if (NormalizePlot) {
  plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[1]]$NormCount)
  plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount / max(plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount)
}

plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2]]$NormCount*-1 #This line flips the values for the second group

#===================
#Output plot
#===================

plot.out <- ggplot(plot.data.bin,aes(x=AAproc,y=NormCount, fill = Group)) + 
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
        axis.text = element_text(size = rel(0.5)),
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
select.samples.binPos <- unique(select.samples.binPos) 
#Due to key, this removes replicates if identical sequence mapped to multiple genes

setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))
select.samples.binPos[,c("BCcount","NormCount","AnimalCount","mainStruct","mismatches"):=
                        list(unlist(lapply(strsplit(paste(BC, collapse=","), ","),function(x) length(unique(x)))),
                             sum(NormCount),
                             unlist(lapply(strsplit(paste(Animals, collapse=","), ","),function(x) length(unique(x)))),
                             paste(unique(structure), collapse=","),
                             median(mismatches)), by=key(select.samples.binPos)]

select.samples.binPos <- unique(select.samples.binPos)
select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount",
                                                  "BCcount","AnimalCount","mainStruct",
                                                  "mismatches"), with = FALSE]

setorder(select.samples.binPos,Group,-NormCount,BCcount,AnimalCount)
setkey(select.samples.binPos,Group)
select.samples.top <- select.samples.binPos[, head(.SD, 20), by=Group]
top.sample <- select.samples.top[J(names(fill.values)[1])]
bottom.sample <- select.samples.top[J(names(fill.values)[2])]
top.sample[,c("Group"):=NULL]
bottom.sample[,c("Group"):=NULL]

out.list <- list(plot=plot.out,
                 plotBin=plot.data.bin,
                 top=top.sample,
                 bottom=bottom.sample)

return(out.list)
}


#'Analyze samples
#'===================

plotPair("totalLib","infectiveLib")$plot

plotPair("CNS100x_Str","totalLib")$plot

out.plot.list <- plotPair("CNS100x_Th","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS100x_Ctx","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS100x_SN","CNS100x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS100x_SN","CNS100x_Th")$plot

plotPair("CNS100x_Ctx","CNS100x_Th")$plot

plotPair("CNS100x_SN","CNS100x_Ctx")$plot

out.plot.list <- plotPair("CNS1000x_Th","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("CNS1000x_Ctx","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

out.plot.list <- plotPair("CNS1000x_SN","CNS1000x_Str")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")

plotPair("CNS1000x_SN","CNS1000x_Th")$plot
plotPair("CNS1000x_Ctx","CNS1000x_Th")$plot
plotPair("CNS1000x_SN","CNS1000x_Ctx")$plot
plotPair("CNS1000x_Str","CNS100x_Str")$plot
plotPair("CNS1000x_Th","CNS100x_Th")$plot
plotPair("CNS1000x_Ctx","CNS100x_Ctx")$plot
plotPair("CNS1000x_SN","CNS100x_SN")$plot

plotPair("CNS100x_Ctx_4wks","CNS100x_Ctx")$plot
plotPair("CNS100x_SN_4wks","CNS100x_SN")$plot
plotPair("CNS100x_Str_4wks","CNS100x_Str")$plot
plotPair("CNS100x_Th_4wks","CNS100x_Th")$plot
plotPair("PerN100x_SC_4wks","PerN100x_Mu_4wks")$plot

# plotPair("CNS1000x_Ctx_4wks","CNS1000x_Ctx")$plot
# plotPair("CNS1000x_SN_4wks","CNS1000x_SN")$plot
# plotPair("CNS1000x_Str_4wks","CNS1000x_Str")$plot
# plotPair("CNS1000x_Th_4wks","CNS1000x_Th")$plot

out.plot.list <- plotPair("PrimN_1000x","PrimN_100x")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

out.plot.list <- plotPair("H293T_1000x","H293T_100x")
out.plot.list$plot
knitr::kable(out.plot.list$top, format = "markdown")
knitr::kable(out.plot.list$bottom, format = "markdown")

plotPair("CNS100x_Str","totalLib",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_Str","CNS100x_Str",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_Th","CNS100x_Th",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_Ctx","CNS100x_Ctx",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_SN","CNS100x_SN",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_SN","CNS1000x_Ctx",size.bin=1,winWidth=4)$plot
plotPair("CNS100x_SN","CNS100x_Ctx",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_Ctx","CNS1000x_Th",size.bin=1,winWidth=4)$plot
plotPair("CNS100x_Ctx","CNS100x_Th",size.bin=1,winWidth=4)$plot
plotPair("CNS1000x_Th","CNS1000x_Str",size.bin=1,winWidth=4)$plot
plotPair("CNS100x_Th","CNS100x_Str",size.bin=1,winWidth=4)$plot
plotPair("PrimN_1000x","PrimN_100x",size.bin=1,winWidth=4)$plot
plotPair("H293T_1000x","H293T_100x",size.bin=1,winWidth=4)$plot

devtools::session_info()
