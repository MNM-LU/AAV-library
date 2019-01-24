# AAV-library
This R-script collection enables identification of fragments and barcodes and allows for barcode quantification
It utilizes knowledge on the surrounding sequence to extract the barcode and requires bbmap, Starcode and Bowtie2 to function

# Download sequencing files
Will be made available at the time of Manuscript publication

# Run command

docker run -d -v /mnt/data/Dropbox\ \(Bjorklund\ Lab\)/Shared/NGS\ data/Sorted\ sequencing\ files/AAVlib-SRA:/home/rstudio/seqFiles -i -p 8787:8787 --name aav-lib bjorklund/aavlib:latest

