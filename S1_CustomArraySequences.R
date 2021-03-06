
#' ---
#' title: "Custom array sequence generation"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' ---

#' This script generates all unique AA sequences for the CustomArray production  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 7.5, fig.height = 8)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)

#'Loading source files
#'===================
source(file.path("functions", "AAtoDNA.R"))
source(file.path("functions", "GeneCodon.R"))
#Override the GeneCodon function with local version containing human codons
unlockBinding("GeneCodon", as.environment("package:GeneGA"))
assign("GeneCodon", GeneCodon, as.environment("package:GeneGA"))

allSequences <- readFasta("input/DNA-lib_RetrogradeTransport.fasta")
AAlist <- data.frame(Class=character(),
                     Family=character(),
                     Strain=character(),
                     Note=character(),
                     Number=character(),
                     Name=character(),
                     AAfragment=character(),
                     stringsAsFactors = FALSE)

#allSequences <- allSequences[124:129] #Debug row


#'Generation of AA table for the selected proteins
#'===================

strt<-Sys.time()
for (i in 1:length(allSequences)){
  thisID <- as.character(ShortRead::id(allSequences[i]))
  thisSeq <- sread(allSequences[i])
  thisAA <- Biostrings::translate(thisSeq, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
  AAlist[i,c("Class","Family","Strain","Note", "Number","Name","AAfragment")] <- c(BBmisc::explode(thisID, sep=","),as.character(thisAA))
}

#'The generateFragments function
#'===================


generateFragments <- function(minLength,maxLength,frequency=1) {
#Generate a table to store all AA sequences
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
    thisFullAA <- AAlist[k,"AAfragment"]
    count = length(fragList[,1])+1
    for (l in seq(1,(width(thisFullAA)-minLength),frequency)){
      for (m in (l+minLength-1):(min(l+maxLength-1,width(thisFullAA)))){
      #Truncate string to the relevant fragment:
      thisFragment <- substr(thisFullAA,l,m) 
      #Take away any sequence that starts with a start codon ATG:
      if (substr(thisFragment,1,1)!="M") { 
      fragList[count,c("Class","Family","Strain","Note","Number",
                       "Name","AAstart","AAstop",
                       "AAfragment")] <- c(AAlist[k,c("Class","Family","Strain",
                                                      "Note","Number","Name")],l,m,thisFragment) 
      ## Inserts the fragment with information into the new data frame 
      count <- count +1
      }
      }
    }
    return(fragList)
  }

fragList <- do.call(rbind,mclapply(1:length(AAlist[,1]),makeAllFrags, 
                                   mc.preschedule = TRUE, mc.cores = detectCores()))

#Control if any sequences contain non-AA characters and save them into a separate list
discardList <- fragList[grep("[[:punct:]|X]",fragList[,"AAfragment"]),]

#Remove any sequence containing non AA characters
fragList <- fragList[grep("[[:punct:]|X]",fragList[,"AAfragment"], invert = TRUE),] 


#Sort the fragments, find unique strings and count number of duplicates
sortedFragments <- rev(sort(table(fragList[,"AAfragment"]))) 

#Run the AAtoDNA function to convert all AA sequences to human codon-optimized DNA sequences
row.names(sortedFragments) <- mclapply(row.names(sortedFragments), fullOPT=FALSE, species="hsa", 
                         AAtoDNA, mc.preschedule = TRUE, mc.set.seed = TRUE,
                         mc.silent = FALSE, mc.cores = detectCores(), mc.cleanup = TRUE) #

sortedFragments <- sortedFragments[order(row.names(sortedFragments))]

return(sortedFragments)
}

#'Execution of the function
#'===================

sortedFragments.14aa <- generateFragments(14,14,1)

#'Add the overhangs for amplication PCR and Gibson assembly into the AAV plasmid
fivePrime <- tolower("AACCTCCAGAGAGGCAACGCT")
threePrime <- tolower("GCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa) <- paste(fivePrime,row.names(sortedFragments.14aa),
                                         threePrime, sep = "")

sortedFragments.14aa.G4S <- generateFragments(14,14,3)
sortedFragments.14aa.A5 <- sortedFragments.14aa.G4S

#'Add the overhangs including G4S spacers for amplication PCR and Gibson assembly into the AAV plasmid

fivePrime <- tolower("AACCTCCAGAGAGGCAACGGAGGCGGAGGAAGT")
threePrime <- tolower("GGAGGCGGCGGAAGCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa.G4S) <- paste(fivePrime,row.names(sortedFragments.14aa.G4S),
                                             threePrime, sep = "")

#'Add the overhangs including A5 spacers for amplication PCR and Gibson assembly into the AAV plasmid

fivePrime <- tolower("AACCTCCAGAGAGGCAACGCTGCTGCAGCAGCC")
threePrime <- tolower("GCAGCTGCAGCTGCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.14aa.A5) <- paste(fivePrime,row.names(sortedFragments.14aa.A5),
                                            threePrime, sep = "")

#'Generate 22aa fragments

sortedFragments.22aa <- generateFragments(22,22,3)

#'Add the overhangs for amplication PCR and Gibson assembly into the AAV plasmid

fivePrime <- tolower("AACCTCCAGAGAGGCAACGCT")
threePrime <- tolower("GCCAGACAAGCAGCTACCGCA")
row.names(sortedFragments.22aa) <- paste(fivePrime,row.names(sortedFragments.22aa),
                                         threePrime, sep = "")

#'Merge all separate fragment lists into one complete list

sortedFragments <- c(sortedFragments.22aa,sortedFragments.14aa,
                     sortedFragments.14aa.A5,sortedFragments.14aa.G4S)

print(paste("Number of unique fragments:",length(unique(names(sortedFragments))), sep=" "))

write.table(c("Sequence",unique(names(sortedFragments))),
            "data/SortedFragments_all.txt",row.names=F,col.names=F,quote=F,sep="\t")

print(Sys.time()-strt)

devtools::session_info()
