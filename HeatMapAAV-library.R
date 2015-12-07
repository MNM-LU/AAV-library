
# HeatMapAAV-library.R creates a color-coded plot over efficiency and specificity in transfection of different parts of assessed proteins. 


# Load all relevant packages, sources all relevant functions. Reads in config-file and imports(?) data file.

library(knitr)
library(ade4)
#library(seqinr)
library(RColorBrewer)
library(ape)
library(gtools)
library(gdata)
library(gplots)
library(ggtree)
library(rncl)
library(ggplot2)
library(geiger)
library(diversitree)
library(dendextend)



# Input specifying what kind of plot to be made? (single gene, avg of genes, ... )

filename <- "~/Dropbox/MNM-Morgan/RNAtablesCompleteBin.rda"
load(filename)

# Define the sites where the data come from. They are ordered the way they eventually will be ordered on the x-axis of the heatmap.

sites <- c("library.table","total.infectiveLib","RatNr1_100x_Str_7","RatNr7_100x_Str_3","RatNr8_100x_Str_11",
           "RatNr15_1000x_Str_15_RatNr19_1000x_Str_22_RatNr20_1000x_Str_24_RatNr21_1000x_Str_19", 
           "RatNr1_100x_Th_8","RatNr7_100x_Th_4","RatNr8_100x_Th_12",
           "RatNr15_1000x_Th_16_RatNr19_1000x_Th_23_RatNr20_1000x_Th_25_RatNr21_1000x_Th_20",
           "RatNr1_100x_SN_5","RatNr7_100x_SN_1","RatNr8_100x_SN_9",
           "RatNr15_1000x_SN_13_RatNr21_1000x_SN_17",
           "RatNr1_100x_Ctx_6","RatNr7_100x_Ctx_2","RatNr8_100x_Ctx_10",
           "RatNr15_1000x_Ctx_14_RatNr19_1000x_Ctx_21_RatNr21_1000x_Ctx_18",
           "primNeuronsNr6_1000x_cDNA_28","primNeuronsNr7_100x_cDNA_29",
           "Cells293Nr2_1000x_cDNA_26","Cells293Nr3_100x_cDNA_27")  # Rename variables eventually to be able to use them as labels in the heatmap.

L <- length(sites)

all.genes <- library.table$V1
all.genes.gsub <- gsub("([,])", "-",gsub("([#])","Number",(gsub("([ ])", "_", all.genes))))

# Loop over number of data tables and put them in a list.

listOfTables <- vector(mode = "list", length = L)
for (i in 1:L){listOfTables[[i]] <- eval(parse(text=sites[i]))} 

# Define a function that merges multiple tables that are in a list.

merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

# Assign the index of the column of the readcounts to the variable "readcount.type". Column 2 contains
# the true readcounts, column 3 the barcode counts and column 4 the normalized readcounts.

readcount.type=4    

# Merge the data tables by gene name, so that missing values don't mess up the order. Replace the NA's with zeros and use the gene names as 
# row names. Finally convert the data table into a matrix so that it can work as input to the "heatmap" function.

merged.df <-  Reduce(function(...) merge(...,by = "V1", all=T), listOfTables)
all.counts.df <- merged.df[ ,seq(from=readcount.type,to=3*L+1,by=3)]
all.counts.df[is.na(all.counts.df)] <- 0
colnames(all.counts.df) <- sites
rownames(all.counts.df) <- all.genes.gsub
all.counts.mat <- data.matrix(all.counts.df) 
# dimnames(all.counts.mat) <- list(all.genes.gsub,sites)

# Rlog-transform



# Import dendrogram, convert to a dendrogram suitable for input to "heatmap".
library(phytools)
tree <- phytools::read.newick(file="tree.phy")


tree <- read_nexus_phylo("treeV2", simplify = FALSE, missing_edge_length = NA) # subtree.phy is a version of tree.phy where "#" has been replaced by "Number".
#tree <- read_newick_phylo("tree.phy", simplify = FALSE, missing_edge_length = NA) # subtree.phy is a version of tree.phy where "#" has been replaced by "Number".
tree.calib <- makeChronosCalib(tree, age.min = 0, age.max = max(tree$edge.length))
chronogram.genes <- chronos(tree) #, calibration = tree.calib


hc.genes <- as.hclust.phylo(chronogram.genes) 
dendrogram.genes <- as.dendrogram(hc.genes)

plot(dendrogram.genes, xlim=c(0,10))


# Do the graphics! label w row dendrogram

# Color scale continuous/binned? "Biased" ?

#colorpalette <- colorRampPalette(heat.colors(n=300,alpha=1),bias=5,alpha=1)
#colorpalette <- heat.colors(n=299,alpha=1)
colorpalette=rainbow(n=300,start=4/6,end=0)

h <- heatmap.2(x = all.counts.mat,
              Colv = NA, 
              Rowv = NA,
            # Rowv = dendrogram.genes,
              dendrogram = "none", #row
            # rowsep = "black",
            # xlab="Brain Region",
            # ylab="Gene fragment",
              key=TRUE,
              keysize=1.5,
              key.title = "Color key",
              density.info = "density",
              tracecol = "slategray2",#3/6,
              key.xlab = attributes(library.table)$names[readcount.type],
              key.ylab = "Count Density",
            # labRow=all.genes,
            # LabCol=sites,
              col = colorpalette)
  