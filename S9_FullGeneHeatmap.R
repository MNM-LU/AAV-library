
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
opts_chunk$set(fig.width = 8, fig.height = 11)
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

#'Setup parameters
#'===================
select.samples <- readRDS("data/normalizedSampleRanges.RDS")
mcols(select.samples)$Sequence <- names(select.samples)
names(select.samples) <- make.names(names(select.samples), unique=TRUE)
length.Table <- data.table(seqnames=names(seqlengths(select.samples)),
                           seqlength=seqlengths(select.samples), key="seqnames")
select.samples <- data.table(as.data.frame(select.samples), key="seqnames")
select.samples[,c("strand","qwidth","cigar","njunc","end"):=NULL]
select.samples <- select.samples[length.Table] #A data.table berge to match seglengths to their respective seqnames

#'Selection of relevant samples
#'===================
select.samples <- select.samples[-grep("4wks|PrimN_1000x_RNA",select.samples$Group),]
select.samples[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(seqnames, ",", fixed=TRUE)]
select.samples[,c("seqnames","Extra"):=NULL]

select.samples.binCat <<- select.samples
setkeyv(select.samples.binCat,c("Group","Category"))
select.samples.binCat[,c("BCcount","NormCount"):=list(unlist(lapply(strsplit(paste(BC, collapse=","), ","),function(x) length(unique(x)))),
                                                      mean(NormCount)), by=key(select.samples.binCat)]
select.samples.binCat <- unique(select.samples.binCat)
select.samples.binCat <- select.samples.binCat[,c("Group","Category","BCcount","NormCount"), with = FALSE]
select.samples.binCat[,totBC:=sum(BCcount), by="Group"]
max.count <- max(select.samples.binCat$totBC)
select.samples.binCat[,BCcountN:=BCcount/totBC*max.count]

length.Table[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(seqnames, ",", fixed=TRUE)]
setkey(length.Table,"Category")
length.Table[,seqlength:=sum(seqlength), by="Category"]
length.Table <- length.Table[,c("Category","seqlength"), with = FALSE]
length.Table <- unique(length.Table)
setkey(select.samples.binCat,"Category")
select.samples.binCat <- select.samples.binCat[length.Table,nomatch=0]
select.samples.binCat[,Category:=gsub("/|_|’","-",Category)]
select.samples.binCat[,BCcountNseq:=BCcountN/seqlength]
select.samples.binCat[,NormCountBC:=BCcountNseq*NormCount]


select.samples.binGene <<- select.samples
setkeyv(select.samples.binGene,c("Group","Category","GeneName"))
select.samples.binGene[,c("BCcount","NormCount"):=list(unlist(lapply(strsplit(paste(BC, collapse=","), ","),function(x) length(unique(x)))),
                                                       mean(NormCount)), by=key(select.samples.binGene)]
select.samples.binGene <- unique(select.samples.binGene)
select.samples.binGene <- select.samples.binGene[,c("Group","GeneName","Category","BCcount","seqlength","NormCount"), with = FALSE]
select.samples.binGene[,GeneName:=gsub("/|_|’","-",GeneName)]
select.samples.binGene[,totBC:=sum(BCcount), by="Group"]
max.count <- max(select.samples.binGene$totBC)
select.samples.binGene[,BCcountN:=BCcount/totBC*max.count]
select.samples.binGene[,BCcountNseq:=BCcountN/seqlength]
select.samples.binGene[,NormCountBC:=BCcountNseq*NormCount]

select.samples.binPos <<- select.samples
setkeyv(select.samples.binPos,c("Group","structure","Sequence"))
select.samples.binPos <- unique(select.samples.binPos) #Due to key, this removes replicates if identical sequence mapped to multiple genes
select.samples.binPos[,AA:=((start+2)/3)+floor((width/3)/2)]
setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))
select.samples.binPos[,c("BCcount","NormCount","AnimalCount","mainStruct","mismatches"):=
                        list(unlist(lapply(strsplit(paste(BC, collapse=","), ","),function(x) length(unique(x)))),
                             mean(NormCount),
                             unlist(lapply(strsplit(paste(Animals, collapse=","), ","),function(x) length(unique(x)))),
                             paste(unique(structure), collapse=","),
                             median(mismatches)), by=key(select.samples.binPos)]
select.samples.binPos <- unique(select.samples.binPos)
select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount","BCcount","AnimalCount","mainStruct","mismatches","seqlength"), with = FALSE]
select.samples.binPos[,totBC:=sum(BCcount), by="Group"]
max.count <- max(select.samples.binPos$totBC)
select.samples.binPos[,BCcountN:=BCcount/totBC*max.count]
select.samples.binPos[,BCcountNanim:=BCcountN*AnimalCount]
select.samples.binPos[,BCcountNseq:=BCcountN/seqlength]
select.samples.binPos[,NormCountBC:=BCcountNanim*NormCount]
select.samples.binPos[,GeneAA:=paste(GeneName," [",AA,"]",sep="")]



#'Plot Heatmaps split by Category
#'===================

plotCategory <- function(select.samples.table,plot.col,sample.select){
setkey(select.samples.table,Group)
select.samples.select <- select.samples.table[sample.select]
eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,")", sep="")))
select.samples.matrix <- acast(select.samples.select, Category~Group, value.var=plot.col) #Utilizes reshape 2 to make matrix for heatmap
select.samples.matrix[is.na(select.samples.matrix)] <- 0
select.samples.matrix <- select.samples.matrix[,sample.select]
return(pheatmap(select.samples.matrix, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE))
}

plotCategory(select.samples.binCat,"BCcountNseq",c("totalLib","CNS100x_Str","CNS100x_Th","CNS100x_Ctx","CNS100x_SN"))
plotCategory(select.samples.binCat,"BCcountNseq",c("CNS1000x_Str","CNS1000x_Th","CNS1000x_Ctx","CNS1000x_SN"))
plotCategory(select.samples.binCat,"BCcountNseq",c("PrimN_100x","PrimN_1000x","293T_100x","293T_1000x"))




#'Plot Heatmaps split by GeneName
#'===================

plotGene <- function(select.samples.table,plot.col,sample.select){
setkey(select.samples.table,Group)
select.samples.select <- select.samples.table[sample.select]
eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,")", sep="")))
select.samples.matrix <- acast(select.samples.select, GeneName~Group, value.var=plot.col) #Utilizes reshape 2 to make matrix for heatmap
select.samples.matrix[is.na(select.samples.matrix)] <- 0
select.samples.matrix <- select.samples.matrix[,sample.select]
return(pheatmap(select.samples.matrix, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE))
}

plotGene(select.samples.binGene,"BCcountNseq",c("totalLib","CNS100x_Str","CNS100x_Th","CNS100x_Ctx","CNS100x_SN"))
plotGene(select.samples.binGene,"BCcountNseq",c("CNS1000x_Str","CNS1000x_Th","CNS1000x_Ctx","CNS1000x_SN"))
plotGene(select.samples.binGene,"BCcountNseq",c("PrimN_100x","PrimN_1000x","293T_100x","293T_1000x"))


#'Selection of top ten fragments per sample
#'===================


plotPos <- function(select.samples.table,plot.col,sample.select){
setkeyv(select.samples.table,"Group")
select.samples.select <- select.samples.table[sample.select]
eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,")", sep="")))
select.samples.topTen <- select.samples.select[, head(.SD, 10), by=Group] 
select.samples.out <- select.samples.select[select.samples.select$GeneAA %in% select.samples.topTen$GeneAA]
select.samples.out <- acast(select.samples.out, GeneAA~Group, value.var=plot.col) #Utilizes reshape 2 to make matrix for heatmap
select.samples.out[is.na(select.samples.out)] <- 0
select.samples.out <- select.samples.out[,sample.select]
pheatmap(select.samples.out, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)
select.samples.topTen[,c("totBC","seqlength","BCcountN","GeneAA","BCcountNseq"):=NULL]
return(knitr::kable(select.samples.topTen, format = "markdown"))
}

plotPos(select.samples.binPos,"BCcountN",c("PrimN_100x","PrimN_1000x","293T_100x","293T_1000x"))
plotPos(select.samples.binPos,"BCcountN",c("totalLib","CNS100x_Str","CNS100x_Th","CNS100x_Ctx","CNS100x_SN"))
plotPos(select.samples.binPos,"BCcountNanim",c("totalLib","CNS100x_Str","CNS100x_Th","CNS100x_Ctx","CNS100x_SN"))
plotPos(select.samples.binPos,"BCcountN",c("CNS1000x_Str","CNS1000x_Th","CNS1000x_Ctx","CNS1000x_SN"))
plotPos(select.samples.binPos,"BCcountNanim",c("CNS1000x_Str","CNS1000x_Th","CNS1000x_Ctx","CNS1000x_SN"))


devtools::session_info()



