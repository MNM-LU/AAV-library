suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(rncl))
suppressPackageStartupMessages(library(data.table))
library(ggplot2)
library(geiger)
library(diversitree)
in.fasta <-  readFasta("DNA Libraries for Retrograde Transport.fasta")
aaSeqs <- Biostrings::translate(sread(in.fasta), genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
name.table <- data.table(as.character(ShortRead::id(in.fasta)),gsub("([ ])", "_",tstrsplit(ShortRead::id(in.fasta), ",", fixed=TRUE)[[6]]))
setnames(name.table, c("V1","V2"), c("FullName", "ShortName"))
names(aaSeqs) <- name.table$ShortName
saveRDS(name.table, file="data/geneNames.rds")
aaLib.file <- tempfile(pattern = "aalib_", tmpdir = tempdir(), fileext = ".fasta")
writeXStringSet(aaSeqs,aaLib.file, format= "fasta")


tree.file <- tempfile(pattern = "results_", tmpdir = tempdir(), fileext = ".tree")

sys.out <- system(paste("~/usearch -cluster_agg ",aaLib.file," -treeout  treeR2.phy -distmxout distance.txt -clusterout clusters.txt -id 0.10 -linkage max "," 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE)
tree <- read_newick_phylo("treeR2.phy", simplify = FALSE, missing_edge_length = NA)
tmp.list <- tree$edge.length
#tmp.list <- as.integer(tmp.list)
tree$edge.length <- tmp.list

tree.calib <- makeChronosCalib(tree, age.min = 0, age.max = max(tmp.list))
dendrogram <- chronos(tree, calibration = tree.calib )

plot(dendrogram)