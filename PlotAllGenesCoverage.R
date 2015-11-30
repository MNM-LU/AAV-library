load("completeLibraryRanges.rda")
suppressPackageStartupMessages(library(ggplot2))
table.analysis <- data.table(as.character(seqnames(complete.ranges)), start(complete.ranges)+(qwidth(complete.ranges)/2), mcols(complete.ranges)$tCount, 1L)
setkey(table.analysis, V1, V2)       
library.table <- table.analysis[,list(ReadCount=sum(V3), BCcount=sum(V4)), by=list(V1, V2)]
library.table$normCount <- library.table$ReadCount/mean(library.table$ReadCount)
library.table[,V2:=(V2+2)/3]
setnames(library.table, c("V1", "V2"), c("GeneName", "AApos"))
min(library.table$AApos)


attach("data/RNAtables.rda")
total.library <- get('RatNr1_100x_Ctx_6_RatNr7_100x_Ctx_2_RatNr8_100x_Ctx_10')

detach("file:data/RNAtables.rda")
setnames(total.library, c("V1", "V2"), c("GeneName", "AApos"))

#library.table$BCcount <- (library.table$BCcount/max(library.table$BCcount))*100
library.table$Library <- "Plasmid"
total.library <- total.library[total.library$ReadCount > 1,]
max(total.library$BCcount)
total.library$BCcount <- (total.library$BCcount/max(total.library$BCcount))*max(library.table$BCcount)
total.library$Library <- "AAV"
total.library <- rbind(total.library,library.table)



ggplot(total.library,aes(x=AApos,y=BCcount))+geom_histogram(size=0.15, binwidth = 10, stat="identity", aes(fill = Library))+theme_bw()+
   
  scale_x_continuous(limit=c(0,3129), breaks=c(seq(0,3000,200)), expand =c(0,0)) +
  scale_fill_manual(name = "Library", values = c("Plasmid" = rgb(157,190,217, maxColorValue = 255),"AAV" = rgb(38,64,135, maxColorValue = 255))) +
  scale_colour_manual(name = "Library", values = c("Plasmid" = rgb(157,190,217, maxColorValue = 255),"AAV" = rgb(38,64,135, maxColorValue = 255))) +
  facet_grid(GeneName~., scales = "free_x", space = "free_x") 

#scale_y_continuous(limit=c(0,550), breaks=c(seq(0,500,250)), expand =c(0,0)) +
# Lighter library alternative color: 193,210,234