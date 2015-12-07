suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))

attach("data/RNAtables.rda")

total.library <- get('library.table')
total.library.2 <- get('RatNr15_1000x_Str_15_RatNr19_1000x_Str_22_RatNr20_1000x_Str_24_RatNr21_1000x_Str_19')
total.library.3 <- get('RatNr1_100x_Str_7_RatNr7_100x_Str_3_RatNr8_100x_Str_11')
detach("file:data/RNAtables.rda")
setnames(total.library, c("V1", "V2"), c("GeneName", "AApos"))
setnames(total.library.2, c("V1", "V2"), c("GeneName", "AApos"))
setnames(total.library.3, c("V1", "V2"), c("GeneName", "AApos"))

#library.table$BCcount <- (library.table$BCcount/max(library.table$BCcount))*100
name.1 <- "Library"
name.2 <- "Str1000x"
name.3 <- "Str100x"


total.library$ReadCount <- log2(total.library$ReadCount)*total.library$BCcount
total.library$ReadCount <- scale(total.library$ReadCount, center = TRUE)
total.library.2$ReadCount <- log2(total.library.2$ReadCount)*total.library.2$BCcount
total.library.2$ReadCount <- scale(total.library.2$ReadCount, center = TRUE)
total.library.3$ReadCount <- log2(total.library.3$ReadCount)*total.library.3$BCcount
total.library.3$ReadCount <- scale(total.library.3$ReadCount, center = TRUE)

o = order(-total.library.2$ReadCount)
total.library.2 <- total.library.2[o]
total.library.2[1:10,]

o = order(-total.library.3$ReadCount)
total.library.3 <- total.library.3[o]
total.library.3[1:10,]

# max(total.library$ReadCount,total.library.2$ReadCount,total.library.3$ReadCount)
# #total.library <- total.library[total.library$ReadCount > 1,]
# total.library$ReadCount <- (total.library$ReadCount/max(total.library$ReadCount))*max(total.library$ReadCount)
# #total.library.2 <- total.library.2[total.library.2$ReadCount > 1,]
# total.library.2$ReadCount <- (total.library.2$ReadCount/max(total.library.2$ReadCount))*max(total.library$ReadCount)
# #total.library.3 <- total.library.3[total.library.3$ReadCount > 1,]
# total.library.3$ReadCount <- (total.library.3$ReadCount/max(total.library.3$ReadCount))*max(total.library$ReadCount)
total.library$Library <- name.1
total.library.2$Library <- name.2
total.library.3$Library <- name.3
total.library <- rbind(total.library, total.library.2, total.library.3) #,total.library.2
total.library[, c("Category", "Protein", "Origin", "Extra", "Number","GeneName") := tstrsplit(GeneName, ",", fixed=TRUE)]
#total.library <- total.library[total.library$GeneName %in% unique(total.library$GeneName)[1:8]]
# total.library$ReadCount <- log2(total.library$ReadCount)
size.bin <- 1
# threshold <- sd(total.library$ReadCount[total.library$Library == name.1])
# listOut <- total.library[total.library$Library != name.1 & total.library$ReadCount > threshold,]
# total.library <- listOut
# nodes <- detectCores()
# cl <- makeCluster(nodes)
# registerDoParallel(cl)
# 
# sapply(total.library, mode)
# #-----------------------
# position <- seq(0,max(total.library$AApos),size.bin)
# total.library <- cbind(total.library, position_bin=cut(total.library$AApos, breaks=position)) 
# total.library <- ddply(total.library, .(position_bin,GeneName,Library,Category,Protein, Origin), summarize,
#                       AApos = position[findInterval(min(AApos),position)],
#                       ReadCount = sum(ReadCount),
#                       BCcount = sum(BCcount), .parallel = TRUE)
# #-----------------------

fill.values <- c("A" = rgb(157,190,217, maxColorValue = 255), "B" = rgb(38,64,135, maxColorValue = 255), "C" = "springgreen4")
names(fill.values) <- c(name.1,name.2,name.3)

ggplot(total.library,aes(x=AApos,y=ReadCount))+geom_histogram(bin=size.bin, stat="identity", aes(fill = Library,y=ReadCount))+theme_bw()+
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  facet_grid(GeneName~., scales = "free_x", space = "free_x") 

#scale_y_continuous(trans=log2_trans()) +
#  scale_x_continuous(limit=c(0,3129), breaks=c(seq(0,3000,200)), expand =c(0,0)) +







#scale_y_continuous(limit=c(0,550), breaks=c(seq(0,500,250)), expand =c(0,0)) +
# Lighter library alternative color: 193,210,234