FROM rocker/rstudio:3.3.2
#install latex for PDF output
RUN apt-get update
RUN apt-get install -y texlive-latex-base texlive-fonts-recommended texlive-latex-extra libcurl4-openssl-dev libxml2-dev lbzip2 libncurses5-dev libncursesw5-dev wget libbz2-dev liblzma-dev curl ncbi-blast+ parallel bowtie2 xvfb openjdk-7-jre-headless openjdk-7-jdk
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools)" -e "install_version('acepack',version = '1.4.1',repos = 'http://cran.us.r-project.org')" \
-e "install_version('ade4', version = '1.7-4', repos = 'http://cran.us.r-project.org')" \
-e "install_version('backports', version = '1.0.4', repos = 'http://cran.us.r-project.org')" \
-e "install_version('BBmisc', version = '1.10', repos = 'http://cran.us.r-project.org')" \
-e "install_version('beanplot', version = '1.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('bitops', version = '1.0-6', repos = 'http://cran.us.r-project.org')" \
-e "install_version('chron', version = '2.3-47', repos = 'http://cran.us.r-project.org')" \
-e "install_version('cluster', version = '2.0.5', repos = 'http://cran.us.r-project.org')" \
-e "install_version('data.table', version = '1.10.0', repos = 'http://cran.us.r-project.org')" \
-e "install_version('doParallel', version = '1.0.10', repos = 'http://cran.us.r-project.org')" \
-e "install_version('evaluate', version = '0.10', repos = 'http://cran.us.r-project.org')" \
-e "install_version('foreign', version = '0.8-67', repos = 'http://cran.us.r-project.org')" \
-e "install_version('formatR', version = '1.4', repos = 'http://cran.us.r-project.org')" \
-e "install_version('Formula', version = '1.2-1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('ggplot2', version = '2.2.0', repos = 'http://cran.us.r-project.org')" \
-e "install_version('gridExtra', version = '2.2.1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('hash', version = '2.2.6', repos = 'http://cran.us.r-project.org')" \
-e "install_version('Hmisc', version = '4.0-0', repos = 'http://cran.us.r-project.org')" \
-e "install_version('htmltools', version = '0.3.5', repos = 'http://cran.us.r-project.org')" \
-e "install_version('hwriter', version = '1.3.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('knitr', version = '1.15.1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('lattice', version = '0.20-34', repos = 'http://cran.us.r-project.org')" \
-e "install_version('latticeExtra', version = '0.6-28', repos = 'http://cran.us.r-project.org')" \
-e "install_version('magrittr', version = '1.5', repos = 'http://cran.us.r-project.org')" \
-e "install_version('Matrix', version = '1.2-7.1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('lattice', version = '0.20-34', repos = 'http://cran.us.r-project.org')" \
-e "install_version('memoise', version = '1.0.0', repos = 'http://cran.us.r-project.org')" \
-e "install_version('munsell', version = '0.4.3', repos = 'http://cran.us.r-project.org')" \
-e "install_version('multicore', version = '0.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('nnet', version = '7.3-12', repos = 'http://cran.us.r-project.org')" \
-e "install_version('packrat', version = '0.4.8-1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('pastecs', version = '1.3-18', repos = 'http://cran.us.r-project.org')" \
-e "install_version('pheatmap', version = '1.0.8', repos = 'http://cran.us.r-project.org')" \
-e "install_version('plyr', version = '1.8.4', repos = 'http://cran.us.r-project.org')" \
-e "install_version('RColorBrewer', version = '1.1-2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('Rcpp', version = '0.12.8', repos = 'http://cran.us.r-project.org')" \
-e "install_version('rmarkdown', version = '1.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('rpart', version = '4.1-10', repos = 'http://cran.us.r-project.org')" \
-e "install_version('scales', version = '0.4.1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('stringdist', version = '0.9.4.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('seqinr', version = '3.3-3', repos = 'http://cran.us.r-project.org')" \
-e "install_version('stringi', version = '1.1.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('stringr', version = '1.1.0', repos = 'http://cran.us.r-project.org')" \
-e "install_version('survival', version = '2.40-1', repos = 'http://cran.us.r-project.org')" \
-e "install_version('VennDiagram', version = '1.6.17', repos = 'http://cran.us.r-project.org')" \
-e "install_version('withr', version = '1.0.2', repos = 'http://cran.us.r-project.org')" \
-e "install_version('yaml', version = '2.1.14', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R')" -e "biocLite()"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R')" -e "biocLite('biovizBase')" \
-e "biocLite('BSgenome')" -e "biocLite('BiocParallel')" -e "biocLite('Biostrings')" \
-e "biocLite('DESeq2')" -e "biocLite('GeneGA')" -e "biocLite('GenomeInfoDb')" \
-e "biocLite('GenomicAlignments')" -e "biocLite('GenomicFeatures')" -e "biocLite('GenomicRanges')" \
-e "biocLite('ggbio')" -e "biocLite('Gviz')" -e "biocLite('IRanges')" \
-e "biocLite('Rsamtools')" -e "biocLite('S4Vectors')" -e "biocLite('ShortRead')" \
-e "biocLite('SummarizedExperiment')" -e "biocLite('XVector')" -e "biocLite('zlibbioc')"
RUN Rscript -e "library(devtools)" -e "install_version('matrixStats', version = '0.51.0', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "library(devtools)" -e "devtools::install_github('guiastrennec/ggplus')"
WORKDIR /root
#install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 -O /root/samtools-1.4.tar.bz2
RUN tar -xf samtools-1.4.tar.bz2
WORKDIR /root/samtools-1.4
RUN make
RUN make install
RUN rm /root/samtools-1.4.tar.bz2
WORKDIR /root
#install Pairfq
RUN curl -sL cpanmin.us | perl - git://github.com/sestaton/Pairfq.git
#install BBmap
WORKDIR /home/rstudio
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_37.02.tar.gz -O BBMmap.tar.gz
RUN tar -xvf BBMmap.tar.gz && rm BBMmap.tar.gz
#install starcode
RUN git clone git://github.com/gui11aume/starcode.git
WORKDIR /home/rstudio/starcode
RUN make
RUN ln -s /home/rstudio/starcode/starcode /usr/bin/starcode
WORKDIR /home/rstudio
# Define working directory.
WORKDIR /fiji
# Install Fiji.
RUN curl -O http://update.imagej.net/bootstrap2.js && jrunscript bootstrap2.js update-force-pristine
# Add fiji to the PATH
ENV PATH $PATH:/fiji
WORKDIR /home/rstudio
#Adding the scripts and environment files
ADD ./AAV-library.Rproj ./commandLineFullRun.sh ./S1_CustomArraySequences.R ./S2_AAV-lib_libraryExtraction.R ./S3_libraryIdentification.R ./S4_Fragment_translation.R ./S5_AAV-lib_tissueExtraction.R ./S6_generateCompleteLibraryRanges.R ./S7_NormalizedGeneIdentification.R ./S8_PlotAllGenesCoverage.R ./S9_FullGeneHeatmap.R ./S10_generateLibAnalysisPlots.R ./S11_slidingMeanTopHits.R ./S12_PlotThreeSampleCoverage.R ./S13_TauPlateAnalysis.R ./
RUN mkdir functions input macros
ADD ./functions/AAtoDNA.R ./functions/GeneCodon.R ./functions/positiveFilter_mRNA.R ./functions/retrieveFASTAQID.R ./functions/
ADD ./input/config.txt ./input/DNA-lib_RetrogradeTransport.fasta ./input/loadlist_excluded.txt ./input/loadlist.txt ./input/Nextera_PCR_for_all_25_tissues_and_5_in_vitro_samples.xlsx ./input/wSet.rda ./input/Dilutions.txt ./input/Grouping.txt ./input/
ADD ./PlaterunnerHeadless.ijm ./macros
RUN mkdir logs data output
RUN chown -R rstudio:rstudio /home/rstudio/*
