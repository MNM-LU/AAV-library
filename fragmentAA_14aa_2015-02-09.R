
#' ---
#' title: "Custom array sequence generation"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This script generates all AA unique AA sequences for the CustomArray production  
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))

#'Loading source files
#'===================
source(file.path("functions", "AAtoDNA.R"))

allSequences <- readFasta("DNA Libraries for Retrograde Transport.fasta")
AAlist <- data.frame(Class=character(),
                     Family=character(),
                     Strain=character(),
                     Note=character(),
                     Number=character(),
                     Name=character(),
                     AAfragment=character(),
                     stringsAsFactors = FALSE)

allSequences <- allSequences[124:129]

strt<-Sys.time()
for (i in 1:length(allSequences)){
  thisID <- as.character(ShortRead::id(allSequences[i]))
  thisSeq <- sread(allSequences[i])
  thisAA <- Biostrings::translate(thisSeq, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
  AAlist[i,c("Class","Family","Strain","Note","Number","Name","AAfragment")] <- c(BBmisc::explode(thisID, sep=","),as.character(thisAA))
}

#'The generateFragments function
#'===================


generateFragments <- function(minLength,maxLength,frequency=1) {
#minLength=14
#maxLength=14

fragList <- data.frame(Class=character(),
                         Family=character(),
                         Strain=character(),
                         Note=character(),
                         Number=character(),
                         Name=character(),
                         AAstart=integer(),
                         AAstop=integer(),
                         AAfragment=character(),
                         stringsAsFactors = FALSE)    
  
  
makeAllFrags <- function(k){
  #k <- 1
    thisFullAA <- AAlist[k,"AAfragment"]
    count = length(fragList[,1])+1
    for (l in seq(1,(width(thisFullAA)-minLength),frequency)){
      for (m in (l+minLength-1):(min(l+maxLength-1,width(thisFullAA)))){
      thisFragment <- substr(thisFullAA,l,m) ##Truncates string to the relevant fragment
      if (substr(thisFragment,1,1)!="M") { ##This takes away any sequence that starts with a start codon ATG
      fragList[count,c("Class","Family","Strain","Note","Number","Name","AAstart","AAstop","AAfragment")] <- c(AAlist[k,c("Class","Family","Strain","Note","Number","Name")],l,m,thisFragment) ## Inserts the fragment with information into the new data frame 
      count <- count +1
      }
      }
    }
    return(fragList)
  }

fragList <- do.call(rbind,mclapply(1:length(AAlist[,1]),makeAllFrags, mc.preschedule = TRUE, mc.cores = detectCores()))

if (length(grep("[[:punct:]|X]",fragList[,"AAfragment"])) > 0) {
discardList <- fragList[grep("[[:punct:]|X]",fragList[,"AAfragment"]),] ##This line saves  any sequence containing non AA characters into a separate list 
fragList <- fragList[grep("[[:punct:]|X]",fragList[,"AAfragment"], invert = TRUE),] ##This line removes any sequence containing non AA characters
}

sortedFragments <- rev(sort(table(fragList[,"AAfragment"]))) ##This line sorts the fragments, finds unique strings and conts number of duplicates

row.names(sortedFragments) <- mclapply(row.names(sortedFragments), fullOPT=FALSE, species="hsa", 
                         AAtoDNA, mc.preschedule = TRUE, mc.set.seed = TRUE,
                         mc.silent = FALSE, mc.cores = detectCores(), mc.cleanup = TRUE) #
return(sortedFragments)
}

#'Execution of the function
#'===================

sortedFragments.14aa <- generateFragments(14,14,1)

fivePrime <- tolower("AACCTCCAGAGAGGCAACGCT")
threePrime <- tolower("GCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa) <- paste(fivePrime,row.names(sortedFragments.14aa),threePrime, sep = "")

sortedFragments.14aa.G4S <- generateFragments(14,14,3)
sortedFragments.14aa.A5 <- sortedFragments.14aa.G4S

fivePrime <- tolower("AACCTCCAGAGAGGCAACGGAGGCGGAGGAAGT")
threePrime <- tolower("GGAGGCGGCGGAAGCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa.G4S) <- paste(fivePrime,row.names(sortedFragments.14aa.G4S),threePrime, sep = "")

fivePrime <- tolower("AACCTCCAGAGAGGCAACGCTGCTGCAGCAGCC")
threePrime <- tolower("GCAGCTGCAGCTGCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa.A5) <- paste(fivePrime,row.names(sortedFragments.14aa.A5),threePrime, sep = "")

sortedFragments.22aa <- generateFragments(22,22,3)
fivePrime <- tolower("AACCTCCAGAGAGGCAACGCT")
threePrime <- tolower("GCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.22aa) <- paste(fivePrime,row.names(sortedFragments.22aa),threePrime, sep = "")
sortedFragments <- c(sortedFragments.14aa,sortedFragments.14aa.G4S,sortedFragments.14aa.A5,sortedFragments.22aa)

write.table(names(sortedFragments),"data/SortedFragments_all_2015-12-25.txt",row.names=F,col.names=F,quote=F,sep="\t")

print(Sys.time()-strt)

devtools::session_info()
