
#' ---
#' title: "Clustering Polypeptide motifs using the Hammock hidden Markov model peptide clustering"
#' author: "Tomas Bjorklund"
#' header-includes: \usepackage{longtable}
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.6in
#' fontsize: 9.5pt
#' ---

#' This script clusters Polypeptide motifs using the Hammock hidden Markov model peptide clustering.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 8, fig.height = 10.5, fig.align = 'center') #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
strt1<-Sys.time()

#'Loading samples
#'===================

all.samples <- readRDS("data/allSamplesDataTable.RDS")
all.samples[,Peptide:= as.character(Peptide),]

setkey(all.samples,Group)

select.samples <- all.samples[J(c("mRNA_30cpc_SN","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_SN","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_30cpc_SN_4wks","mRNA_30cpc_Th_4wks","mRNA_30cpc_Ctx_4wks","mRNA_3cpc_SN_4wks","mRNA_3cpc_Th_4wks","mRNA_3cpc_Ctx_4wks"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Animalcount:=as.integer(mclapply(Animals, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount+Animalcount-1,]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/trspSamplesPeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)



#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/Hammock", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/trspSamplesPeptides.fasta -d /home/rstudio/data/Hammock -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
hammock.log <- data.table(readLines("data/Hammock/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)


#'Generation of Weblogo visualization
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/Hammock/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))
id.order <- as.list(ham.clusters$cluster_id)
ham.clusters.all <- data.table(read.table("/home/rstudio/data/Hammock/final_clusters_sequences.tsv", header = TRUE, skip = 0, sep="\t",
                                          stringsAsFactors = FALSE, fill=TRUE))
setkey(select.samples,Peptide)
setkey(select.samples.trsp,Peptide)

unlink("/home/rstudio/data/WEBlogos", recursive = TRUE, force = FALSE)
dir.create(file.path("/home/rstudio/data/", "WEBlogos"), showWarnings = FALSE)

setkey(ham.clusters.all,cluster_id)
setkey(ham.clusters,cluster_id)




generateWeblogo <- function(in.name) {
  #in.name <- ham.clusters$cluster_id[2]
  tmp <- system(paste("weblogo --format PDF --sequence-type protein --size large --errorbars NO --resolution 299 --composition equiprobable --color-scheme chemistry < /home/rstudio/data/Hammock/alignments_final/", in.name, ".aln > /home/rstudio/data/WEBlogos/",in.name,".pdf", sep = ""),
                intern = TRUE, ignore.stdout = FALSE)
  this.main <- ham.clusters[J(in.name)]
  main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]
  cat('\n')
  cat("## Peptide",this.main$main_sequence,"from",main.gene,"with cluster number",in.name, "\n")
  cat('\n')
  cat('\n')
  cat(paste0("![Peptide: ",this.main$main_sequence," from ",main.gene," with cluster number ",in.name,"](/home/rstudio/data/WEBlogos/",in.name,".pdf)"))
  cat('\n')
  print(knitr::kable(this.main[,c(1:8,10)], format = "latex") %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  cat('\n')
  this.cluster <- ham.clusters.all[J(in.name)]
  out.table <- knitr::kable(this.cluster[,c(1:7)], format = "latex")
  print(column_spec(out.table, 3, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  cat('\n')
  this.found <- select.samples[J(this.cluster$sequence)]
  cat('\n')
  out.table <- knitr::kable(this.found[,c(20,15,3,5,4,1,23)], format = "latex")
  print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  cat('\n')
  cat("\n\n\\pagebreak\n")
  cat("\n\n\\clearpage\n")
}

#+ results = 'asis'
lapply(id.order, generateWeblogo)
#+ results = 'markup'


# setkey(select.samples.trsp,Peptide)
# select.samples.trsp.select <- select.samples.trsp[J(c("PPDELNLTTASLPL"))]



print("Total analysis time:")
print(Sys.time()-strt1)

devtools::session_info()

