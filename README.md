# Instructions on how to run the analysis pipeline
This repository contains the complete analysis pipeline from the paper "A novel systematic capsid evolution approach performed in vivo for the
design of AAV vectors with tailored properties and tropism" by Davidssson et al., published in PNAS 2019 (in press).

The R-script collection enables the identification of fragments and barcodes and allows for barcode quantification.
It utilizes knowledge on the surrounding sequence to extract the barcode and requires bbmap, Starcode, and Bowtie2 to function.
The depository is structured around a Dockerfile, which enables the build of a self-contained encapsulation of all required dependencies, including R-studio, bowtie2 et.c., The version presented in the paper can be downloaded pre-built from Docker Hub (See below).

# Sequencing datasets
The raw sequencing files have been uploaded to the NCBI SRA with the Accession number PRJNA473475 ( https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473475 ). Unfortunately, there is no simple command to download all sequencing files from there at once and retaining the original filenames. Thus, to replicate the findings of the paper, it is recommended that you download the entire sequencing dataset as a single (9.4GB) compressed file using the following command (in Linux/MacOS):

```
wget -qO- https://www.dropbox.com/s/6wugrf85ekffqdc/AAVlib-SRA.tar.gz | tar -xvz
```

This generates a folder in the current working directory called "AAVlib-SRA"

If downloaded through other means or placed in a different directory, please modify this path in the commands below accordingly.

It can be viewed online as a directory of files here:
https://www.dropbox.com/sh/ijwbt8nj28hws2q/AADwrXqehtoxxbq-zRLXW6-Da?dl=0

# Pre-built Docker image
To ease execution, the entire workflow has been encapsulated in a Docker container with a web-based interface through Rstudio. This container is publicly available on Docker Hub as bjorklund/aavlib:v0.2 (https://hub.docker.com/r/bjorklund/aavlib )
To execute it on any computer with a functioning docker installed, execute the following string:
```
docker run -d -v "$(pwd)"/AAVlib-SRA:/home/rstudio/seqFiles -i -p 9797:8787 --name aav-lib bjorklund/aavlib:v0.2
```

The string “/path-to_sequenceFolder” needs to be updated with the full path to the folder housing the downloaded sequencing files above.

# Building the Docker container from the git repository

Clone the repository using git:
```
git clone https://bitbucket.org/MNM-LU/aav-library.git
```
Enter the generated directory using:
```
cd aav-library
```

## Adding personalized usearch license

The usearch package, required for the generation of the phylotree in script S10 is a personal non-transferable license. Thus, to build in into the docker container, visit https://www.drive5.com/usearch/download.html and select to download version 10.0.240 for Linux. (Later versions could also work but are not tested.)
In an email, you will receive a URL with a personal license number. e.g., https://drive5.com/cgi-bin/upload3.py?license=20191113089230462 (the number is fake here). Copy the number and paste it into the Dockerfile in this repository where it is stated "MY_LICENSE_HERE". After this, it should build with all functions.  

Build the docker container using
```
docker build -t aavlib:local .
```
Please allow for at least an hour of compilation time, depending on the hardware used.
Then leave the container directory:

```
cd ..
```

Execusion of the docker is then almost identical to the pre built container:

```
docker run -d -v "$(pwd)"/AAVlib-SRA:/home/rstudio/seqFiles -i -p 9797:8787 --name aav-lib aavlib:local
```

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