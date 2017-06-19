
#' ---
#' title: "Tau platerunner analysis"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This script coregisters the Platerunner samples from Tau competition assay using ImageJ and plots the samples.  
suppressPackageStartupMessages(library(knitr))
opts_chunk$set(fig.width = 8, fig.height = 6.5)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))

#Analysis of Tau Platerunner data
system2("xvfb-run", args = "-a /fiji/ImageJ-linux64 -macro /home/rstudio/macros/PlaterunnerHeadless.ijm", stdout = "logs/log_TauPlaterunner.txt")


GFPdata <- read.table("output/NeuronTau_inVitro_GFP_Results.xls", header = TRUE, skip = 0, sep="\t")
mCherrydata <- read.table("output/NeuronTau_inVitro_mCherry_Results.xls", header = TRUE, skip = 0, sep="\t")
mCherrydata$Area <- mCherrydata$Min <- mCherrydata$Max <- GFPdata$Area <- GFPdata$Max <- GFPdata$Min <- NULL
totalData <- merge(GFPdata,mCherrydata,by=c("X"))
colnames(totalData) <- c("Image","GFP","mCherry")
totalData$GMratio <- totalData$GFP/totalData$mCherry

Grouping <- read.table("input/Grouping.txt", header = TRUE, sep="\t")
Dilutions <- read.table("input/Dilutions.txt", header = TRUE, sep="\t")

Plates <- seq(0,960,96)
totalData <- cbind(totalData, Plate=cut(totalData$Image, breaks=Plates)) 
levels(totalData$Plate) <- 1:10
totalData$Well <- totalData$Image - (96*(as.integer(totalData$Plate)-1))
Row <- seq(0,96,12)
totalData <- cbind(totalData, Row=cut(totalData$Well, breaks=Row)) 
levels(totalData$Row) <- 1:12
totalData$Column <- totalData$Well - (12*(as.integer(totalData$Row)-1))
totalData <- merge(totalData,Grouping,by=c("Plate"))
Virus <- seq(0,8,4)
totalData <- cbind(totalData, Virus=cut(as.integer(totalData$Row), breaks=Virus)) 
levels(totalData$Virus) <- c("MNM009","MNM017")
totalData <- merge(totalData,Dilutions,by=c("Column"))
totalData <- totalData[totalData[,"Image"] != 510 & totalData[,"Image"] != 824,]



totalData <- data.table(totalData)
totalData[totalData[,.I[Concentration==0]],"GFPblank":=mean(GFP),by=c("Species","Virus","Concentration")]
totalData[,"GFPblank":=max(GFPblank, na.rm = TRUE),by=c("Species","Virus")]
totalData[totalData[,.I[Concentration==0]],"mCherryblank":=mean(mCherry),by=c("Species","Virus","Concentration")]
totalData[,"mCherryblank":=max(mCherryblank, na.rm = TRUE),by=c("Species","Virus")]
totalData[,"GFPtoBlank":=(GFP/GFPblank)]
totalData[,"mCherrytoBlank":=(mCherry/mCherryblank)]
totalData[,"GMratioToBlank":=(GFPtoBlank/mCherrytoBlank)]
totalDataSums <- totalData
totalDataSums[,c("GFPmean", 
             "GFP_SEM", 
             "mCherrymean", 
             "mCherry_SEM", 
             "GMratioMean", 
             "GMratio_SEM",
             "GFPtoBlankMean",
             "GFPtoBlankSEM",
             "mCherrytoBlankMean",
             "mCherrytoBlankSEM",
             "GMratioToBlankMean",
             "GMratioToBlankSEM"):=list(mean(GFP),
                                         sd(GFP)/sqrt(length(GFP)),
                                         mean(mCherry),
                                         sd(mCherry)/sqrt(length(mCherry)),
                                         mean(GMratio),
                                         sd(GMratio)/sqrt(length(GMratio)),
                                         mean(GFPtoBlank),
                                         sd(GFPtoBlank)/sqrt(length(GFPtoBlank)),
                                         mean(mCherrytoBlank),
                                         sd(mCherrytoBlank)/sqrt(length(mCherrytoBlank)),
                                         mean(GMratioToBlank),
                                         sd(GMratioToBlank)/sqrt(length(GMratioToBlank))
                                         ), by=c("Species","Virus","Concentration")]


setorder(totalDataSums,Image)
totalDataSums <- unique(totalDataSums, by=c("Species","Virus","Concentration"))
totalDataSums[, c("GFP","mCherry", "GMratio", "GFPblank", "mCherryblank", "GFPtoBlank", "mCherrytoBlank","GMratioToBlank"):=NULL]

plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM009"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<12]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM009 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))

plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM009"]]]
plotData <- plotData[plotData[,.I[Species=="T39" | Species=="BSA" | Species=="K18"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<12]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM009 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))

plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM009"]]]
plotData <- plotData[plotData[,.I[Species=="T39" | Species=="T44" | Species=="K18"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<12]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM009 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))


plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM017"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<11]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM017 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))

plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM017"]]]
plotData <- plotData[plotData[,.I[Species=="T39" | Species=="BSA" | Species=="K18"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<12]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM017 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))


plotData <- totalDataSums[totalDataSums[,.I[Virus=="MNM017"]]]
plotData <- plotData[plotData[,.I[Species=="T39" | Species=="T44" | Species=="K18"]]]
plotData <- plotData[plotData[,.I[Column>1 & Column<11]]]

d <- ggplot(plotData, aes(x=Concentration, y=GMratioMean, group=Species))

d + geom_ribbon(aes(ymin=GMratioMean-GMratio_SEM, ymax=GMratioMean+GMratio_SEM),  alpha=0.08) + 
  geom_line(aes(colour = factor(Species))) + labs(title = "MNM017 to MNM025 ratio", y = "GFP to mCherry ratio", x = "Tau concentration") + 
  theme_bw(base_size = 12, base_family = "") +
  theme(axis.line = element_line(size = 0.4))  + coord_trans(x = "log10") +
  scale_x_continuous(breaks=c(1e-10,1e-9,1e-8,1e-7, 1e-6,1e-5), expand =c(0,0))



