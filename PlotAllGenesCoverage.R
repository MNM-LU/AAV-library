suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(lumi))

load("completeLibraryRanges.rda")
table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L)
setkey(table.analysis, V1, V2)       
library.table <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4), ReadSD=sd(V3)), by=list(V1, V2)]
library.table$normCount <- library.table$ReadCount/mean(library.table$ReadCount)
library.table[,V2:=(V2+2)/3]
setnames(library.table, c("V1", "V2"), c("GeneName", "AApos"))
min(library.table$AApos)


attach("data/RNAtables.rda")

total.library <- get('total.infectiveLib')
total.library.2 <- get('RatNr15_1000x_Str_15_RatNr19_1000x_Str_22_RatNr20_1000x_Str_24_RatNr21_1000x_Str_19')
detach("file:data/RNAtables.rda")
setnames(total.library, c("V1", "V2"), c("GeneName", "AApos"))
setnames(total.library.2, c("V1", "V2"), c("GeneName", "AApos"))
corr.table <- merge(library.table, total.library, by= c("GeneName", "AApos"))


#library.table$BCcount <- (library.table$BCcount/max(library.table$BCcount))*100
name.1 <- "Plasmid"
name.2 <- "AAV"
name.3 <- ""
library.table$Library <- name.1

max(total.library$BCcount)
#total.library <- total.library[total.library$ReadCount > 1,]
total.library$BCcount <- (total.library$BCcount/max(total.library$BCcount))*max(library.table$BCcount)
#total.library.2 <- total.library.2[total.library$ReadCount > 1,]
total.library.2$BCcount <- (total.library.2$BCcount/max(total.library.2$BCcount))*max(library.table$BCcount)
total.library$Library <- name.2
total.library.2$Library <- name.3
total.library <- rbind(library.table,total.library) #,total.library.2
total.library[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
corr.table[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
#total.library <- total.library[total.library$GeneName %in% unique(total.library$GeneName)[1:8]]

nodes <- detectCores()
cl <- makeCluster(nodes)
registerDoParallel(cl)

size.bin <- 1
sapply(total.library, mode)
#-----------------------
position <- seq(0,max(total.library$AApos),size.bin)
total.library <- cbind(total.library, position_bin=cut(total.library$AApos, breaks=position)) 
total.library <- ddply(total.library, .(position_bin,GeneName,Library,Category,Protein, Origin), summarize,
                      AApos = position[findInterval(min(AApos),position)],
                      ReadCount = sum(ReadCount),
                      BCcount = sum(BCcount), .parallel = TRUE)
#-----------------------

fill.values <- c("A" = rgb(157,190,217, maxColorValue = 255), "B" = rgb(38,64,135, maxColorValue = 255), "C" = "springgreen4")
names(fill.values) <- c(name.1,name.2,name.3)

ggplot(total.library,aes(x=AApos,y=BCcount))+geom_histogram(bin=size.bin, stat="identity", aes(fill = Library,y=BCcount))+theme_bw()+
  scale_x_continuous(limit=c(0,3129), breaks=c(seq(0,3000,200)), expand =c(0,0)) +
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  facet_grid(GeneName~., scales = "free_x", space = "free_x") 

corr.table[, c("vstCount.y") := vst(ReadCount.y, ReadSD.y)]
p <- ggplot(data = corr.table, aes(x = vstCount.x, y = vstCount.y)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() +
  scale_x_continuous(trans="log2") + 
  scale_y_continuous(trans="log2") 

, expand =c(0,0)
qplot(corr.table[,c("vstCount.x",vstCount.y)])


vst(u, std






#scale_y_continuous(limit=c(0,550), breaks=c(seq(0,500,250)), expand =c(0,0)) +
# Lighter library alternative color: 193,210,234