seqs.original <- readFasta("DNA Libraries for Retrograde Transport.fasta")

seqs.AA <- Biostrings::translate(sread(seqs.original), genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

source("AAtoDNA.R")
seqs.optimized = ShortRead(DNAStringSet(sapply(seqs.AA, function(x) AAtoDNA(x, species="hsa"))), BStringSet(gsub("([ ])", "_", ShortRead::id(seqs.original))))

writeFasta(seqs.optimized,"libraryIndex.fa")
