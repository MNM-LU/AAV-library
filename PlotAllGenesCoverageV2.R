suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))



all.samples <- readRDS("data/normalizedSampleRanges.RDS")
# fill.values <- c("CNS1000x_Str" = rgb(157,190,217, maxColorValue = 255), 
#                  "CNS1000x_Ctx" = rgb(38,64,135, maxColorValue = 255), 
#                  "CNS1000x_SN" = "springgreen4")

fill.values <- c("totalLib" = rgb(38,64,135, maxColorValue = 255), 
                 "infectiveLib_V2" = rgb(157,190,217, maxColorValue = 255))

total.AAV.samples <- all.samples[!(mcols(all.samples)$Group %in% c("totalLib","infectiveLib"))]
mcols(total.AAV.samples)$Group <- "infectiveLib"
all.samples <- append(all.samples,total.AAV.samples)
select.samples <- all.samples[mcols(all.samples)$Group %in% names(fill.values)]
trim.names <- data.table(GeneName=levels(seqnames(select.samples)))
trim.names[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
levels(seqnames(select.samples)) <- trim.names$GeneName



#total.library$ReadCount <- scale(total.library$ReadCount, center = TRUE)
# o = order(-total.library.2$ReadCount)
# total.library.2 <- total.library.2[o]
# total.library.2[1:10,]
# threshold <- sd(total.library$ReadCount[total.library$Library == name.1])
# listOut <- total.library[total.library$Library != name.1 & total.library$ReadCount > threshold,]
# total.library <- listOut


geneTable <- data.frame(seqlengths(select.samples))
colnames(geneTable) <- "SeqLength"
geneTable$GeneName <- row.names(geneTable)

plot.data <- data.frame(GeneName=seqnames(select.samples),
                        Group=factor(mcols(select.samples)$Group, levels = names(fill.values)),
                        AApos=start(select.samples)+(width(select.samples)/2),
                        BC=mcols(select.samples)$BC,
                        RNAcount=mcols(select.samples)$RNAcount,
                        NormCount=mcols(select.samples)$NormCount)
plot.data <- merge(plot.data,geneTable)
plot.data$AAproc <- (plot.data$AApos/plot.data$SeqLength)*100

levels(plot.data$Group)

# nodes <- detectCores()
# cl <- makeCluster(nodes)
# registerDoParallel(cl)

size.bin <- 2
#-----------------------
position <- seq(0,100,size.bin)
plot.data <- cbind(plot.data, position_bin=cut(plot.data$AAproc, breaks=position)) 
plot.data.bin <- ddply(plot.data, .(position_bin,GeneName,Group), summarize,
                      AAproc = position[findInterval(min(AAproc),position)],
                      NormCount = log2(mean(RNAcount)+1)*mean(BC)) # 
#-----------------------                           
plot.data.bin[plot.data.bin$Group == names(fill.values)[2],"NormCount"] <- plot.data.bin[plot.data.bin$Group == names(fill.values)[2],"NormCount"]*-1

ggplot(plot.data.bin,aes(x=AAproc,y=NormCount, fill = Group))+geom_bar(stat="identity", position="identity")+theme_bw()+
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  scale_x_continuous(limit=c(0,100), breaks=c(seq(0,100,20)), expand =c(0,0)) +
  facet_wrap(~ GeneName, ncol=5)
  #facet_grid(GeneName~., scales = "free_x", space = "free_x") 

#scale_y_continuous(trans=log2_trans()) +
#  scale_x_continuous(limit=c(0,3129), breaks=c(seq(0,3000,200)), expand =c(0,0)) +







#scale_y_continuous(limit=c(0,550), breaks=c(seq(0,500,250)), expand =c(0,0)) +
# Lighter library alternative color: 193,210,234