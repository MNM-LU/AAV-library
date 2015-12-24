require(ShortRead)
require(foreach)
require(parallel)
require(GeneGA)
setwd("~/Dropbox (Bjorklund Lab)/mnm group files/AAV WGA project/R analysis")
source("AAtoDNA.R")
minLength <- 14
maxLength <- 14
allSequences <- readFasta("DNA Libraries for Retrograde Transport.fasta")
#i=1
AAlist <- data.frame(NA,NA,NA,NA,NA,NA,NA)
colnames(AAlist) <- c("Class","Family","Strain","Note","Number","Name","AAfragment")
fragList <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
colnames(fragList) <- c("Class","Family","Strain","Note","Number","Name","AAstart","AAstop","AAfragment")
strt<-Sys.time()
for (i in 1:length(allSequences)){
#for (i in 1:5){
  #allSequences[i]
  thisID <- as.character(id(allSequences[i]))
  thisSeq <- sread(allSequences[i])
  thisAA <- Biostrings::translate(thisSeq, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
  AAlist[i,c("Class","Family","Strain","Note","Number","Name","AAfragment")] <- c(BBmisc::explode(thisID, sep=","),as.character(thisAA))
}

for (k in 1:length(AAlist[,1])){
  thisFullAA <- AAlist[k,"AAfragment"]
  count = length(fragList[,1])+1
  for (l in 1:(width(thisFullAA)-minLength)){
    for (m in (l+minLength):(min(l+maxLength,width(thisFullAA)))){
    thisFragment <- substr(thisFullAA,l,m) ##Truncates string to the relevant fragment
    if (substr(thisFragment,1,1)!="M") { ##This takes away any sequence that starts with a start codon ATG
    fragList[count,c("Class","Family","Strain","Note","Number","Name","AAstart","AAstop","AAfragment")] <- c(AAlist[k,c("Class","Family","Strain","Note","Number","Name")],l,m,thisFragment) ## Inserts the fragment with information into the new data frame 
    count <- count +1
    }
    }
  }
}
fragList <- fragList[2:length(fragList[,1]),] #removes the NA line
discardList <- fragList[grep("[[:punct:]]",fragList[,"AAfragment"]),] ##This line saves  any sequence containing non AA characters into a separate list 
fragList <- fragList[grep("[[:punct:]]",fragList[,"AAfragment"], invert = TRUE),] ##This line removes any sequence containing non AA characters
fivePrime <- tolower("TTATGGTACCAtgggT")
threePrime <- tolower("agtggTGGCGGAGGTT")
selected <- fragList[,"AAfragment"]
sortedFragments <- rev(sort(table(selected))) ##This line sorts the fragments, finds unique strings and conts number of duplicates

sequenceList <- mclapply(row.names(sortedFragments), fullOPT=FALSE, species="hsa", AAtoDNA, mc.preschedule = TRUE, mc.set.seed = TRUE,
mc.silent = FALSE, mc.cores = (detectCores()), mc.cleanup = TRUE)

row.names(sortedFragments) <- paste(fivePrime,sequenceList,threePrime, sep = "")

# library(GeneGA)
# source("AAtoDNA.R")
#tempAA <- "MEPYTDSNVPEGDSVPEAHWTKLQHSLDTALRRARSAPAEPIAARVTG"
#tempDNA <- AAtoDNA(tempAA,fullOPT=TRUE, species="hsa", optIt=1)
#tempDNA2 <- AAtoDNA(tempAA,fullOPT=FALSE, species="hsa")
#tempDNA == tempDNA2
# backAA <- Biostrings::translate(DNAString(tempDNA), genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
save(AAlist, file = "GeneList.RData")
save(fragList, file = "FragmentList.RData")
save(discardList, file = "DiscardedList.RData")
write.table(sortedFragments,"SortedFragments.txt",row.names=T,col.names=F,quote=F,sep="\t")
print(Sys.time()-strt)
