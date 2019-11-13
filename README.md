# Instructions on how to run the analysis pipeline
This repository contains the complete analysis pipeline from the paper "A novel systematic capsid evolution approach performed in vivo for the
design of AAV vectors with tailored properties and tropism" by Davidssson et al., published in PNAS 2019 (in press9.

The R-script collection enables identification of fragments and barcodes and allows for barcode quantification
It utilizes knowledge on the surrounding sequence to extract the barcode and requires bbmap, Starcode and Bowtie2 to function.
The depository is structured around a Dockerfile which enables the build of a self-contained encapulation of all required dependencies, including R-studio, bowtie2 et.c., The version presented in the paper can be downloaded pre-built from Docker Hub (See below).

# Sequencing datasets
The raw sequencing files have been uploaded to the NCBI SRA with the Accession number PRJNA473475 ( https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473475 ). These data will be publicly released at the point of publication. Until then, the entire sequencing dataset can be downloaded as a single (9.4GB) compressed file using the following command (in Linux/MacOS):

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
To ease execution, the entire workflow has been encapsulated in a Docker container with a web-based interface through Rstudio. This container is publicly available on Docker Hub as bjorklund/aavlib:v0.2 (https://hub.docker.com/r/bjorklund/aavlib )

To execute it on any computer with a functioning docker installed, execute the following string:
```
docker run -d -v /path-to_sequenceFolder:/home/rstudio/seqFiles -i -p 9797:8787 --name aav-lib bjorklund/aavlib:v0.2
```
The string “/path-to_sequenceFolder” needs to be updated with the full path to the folder housing the downloaded sequencing files above. 

# Running the workflow
Once this Docker container is running, an instance of Rstudo is available on port 9797 (This is an arbitrary port and can be changed to any free port in the docker command above). If executed on the same computer, then point any modern web browser to http://localhost:9797
The username:password for Rstudio is rstudio:rstudio

To execute all workflows and generate a PDF for each dataset, execute all lines in the R-script “S0_RunAll.R”. This will render PDF files in the “output” folder. To execute the code line-by-line, you will need to run the R-scripts S1 to S15 in this folder in the correct order as they generate dependencies for each other.

Note: The complete output can be downloaded from the following link and is called “Davidsson et al., S14 bioinformatics.pdf”.
https://www.dropbox.com/s/gioe1l3cju04lpa/Davidsson%20et%20al.%2C%20S14%20bioinformatics.pdf?dl=1
This will also be part of the online supplement at PNAS.

# Recommended hardware
The original analysis was performed on a “fat node” Linux server with dual Intel Xenon E5-2650 v2 (16 cores) and 256 Gb RAM running Debian 9 and docker. On this system the complete analysis takes approx. 12 hours. We have successfully run the complete analysis on an 8-core server and 64 Gb RAM but we cannot guarantee successful execution on hardware with less memory. 

# Source code and Complete data output

The datasets required to re-run this analysis pipeline are available in the NCBI Sequence Read Archive (SRA) with the accession number PRJNA473475 ( https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473475 ). The R-based workflow is publicly available as a Git repository at https://bitbucket.org/MNM-LU/aav-library with the dockerfile to generate a docker image. A pre-built Docker image is also available on Docker Hub named: Bjorklund/aavlib:v0.2. The output from the entire bioinformatics pipeline is available at https://www.dropbox.com/s/9p14lyjs6e5gt9w/Bioinformatics-output.pdf?dl=1 contains the complete formatted output of the bioinformatics pipeline. This will be made publicly available as part of the Git repository at the time of publication.

