AAtoDNA <- function(inAA,species="ec",fullOPT=FALSE,optIt=1){
  require(Hmisc)
  library(GeneGA)
  humanCodon <- read.table(header = TRUE, 
                           stringsAsFactors = FALSE, 
text="AA DNA
A   gcc
C	tgc
D	gac
E	gag
F	ttc
G	ggc
H	cac
I	atc
K	aag
L	ctg
M	atg
N	aac
P	ccc
Q	cag
R	cgg
S	agc
T	acc
V	gtg
W	tgg
Y	tac")

inAA <- toupper(as.character(inAA))

outDNA <- toupper(sedit(inAA,humanCodon$AA, humanCodon$DNA, wild.literal=FALSE))
outDNA <- DNAString(outDNA)
outDNA <- GeneCodon(as.character(outDNA) , organism = species)
#wSet <- read.table("~/Dropbox (Bjorklund Lab)/mnm group files/AAV WGA project/R analysis/OrganismTable.txt", row.names = 1, header = TRUE, skip = 0, sep="\t")

if (fullOPT == TRUE){
   GeneGAoutput <- GeneGA(sequence = outDNA, popSize = 50, iters =optIt, crossoverRate = 0.2,
                   mutationChance = 0.05, region = NULL, organism = species, 
                   showGeneration = TRUE, frontSeq = NULL, ramp=FALSE, numcode=1)
   outDNA <- GeneGAoutput
}


  outDNA
}
