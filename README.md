# Instructions on how to run the analysis pipeline
This repository contains the complete analysis pipeline from the paper "A novel systematic capsid evolution approach performed in vivo for the
design of AAV vectors with tailored properties and tropism" by Davidssson et al., published in PNAS 2019 (in press9.

The R-script collection enables identification of fragments and barcodes and allows for barcode quantification
It utilizes knowledge on the surrounding sequence to extract the barcode and requires bbmap, Starcode and Bowtie2 to function.
The depository is structured around a Dockerfile which enables the build of a self-contained encapulation of all required dependencies, including R-studio, bowtie2 et.c., The version presented in the paper can be downloaded pre-built from Docker Hub (See below).

# Sequencing datasets
The raw sequencing files have been uploaded to the NCBI SRA with the project number SRP149133. These data will be publicly released at the point of publication. Until then, the entire sequencing dataset can be downloaded as a single (9.4GB) compressed file using the following command (in Linux/MacOS):

```
wget https://www.dropbox.com/s/6wugrf85ekffqdc/AAVlib-SRA.tar.gz 
```

It can be viewed online as a directory of files here:
https://www.dropbox.com/sh/ijwbt8nj28hws2q/AADwrXqehtoxxbq-zRLXW6-Da?dl=0

This file needs to be extracted into a folder somewhere on the computer where to the path is referred to below as “/path-to_sequenceFolder”.

# Pre-built Docker image
To ease execution, the entire workflow has been encapsulated in a Docker container with a web-based interface through Rstudio. This container is publicly available on Docker Hub as bjorklund/aavlib:v0.2 (https://hub.docker.com/r/bjorklund/aavlib )
To execute it on any computer with a functioning docker installed, execute the following string:
```
docker run -d -v /path-to_sequenceFolder:/home/rstudio/seqFiles -i -p 9797:8787 --name aav-lib bjorklund/aavlib:v0.2
```

The string “/path-to_sequenceFolder” needs to be updated with the full path to the folder housing the downloaded sequencing files above.

# Building the Docker container from the git repository

