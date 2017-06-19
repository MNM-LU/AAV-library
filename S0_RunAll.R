# Generation of a lookup table using the Pacific Biosciences RSII CCS reads of full length PCR free plasmid fragmens
sys.args <- paste("-e",shQuote("rmarkdown::render('DNA_LibMapping.R')"), "$PWD config_Syn9-10_Lib_PacBio.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_Syn9-10_Lib_PacBio.txt")
system2("mv", args="DNA_LibMapping.pdf output/DNA_Syn9-10_Lib_PacBio.pdf")

#Mapping of stable clones from sincle cell sorting to the PacBio RSII LUT and map with RNA-editing efficacy
sys.args <- paste("-e",shQuote("rmarkdown::render('mRNA-spliceCount_SingleCell.R')"), "$PWD config_singleCell_Syn9-10_stable.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_singleCell_Syn9-10_stable.txt")
system2("mv", args="mRNA-spliceCount_SingleCell.pdf output/mRNA-spliceCount_singleCell_Syn9-10_stable.pdf")

# Generation of LUT from from both PacBio and emulsion PCR from SalI digested and ligated plasmid library. With the use of a larger sequencing dataset, more barcodes can be recovered.
sys.args <- paste("-e",shQuote("rmarkdown::render('DNA_LibMapping_V2.R')"), "$PWD config_pSyn9-10_totalLib.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_pSyn9-10_totalLib.txt")
system2("mv", args="DNA_LibMapping_V2.pdf output/pSyn9-10_totalLib.pdf")

# Evaluation of RNA-editing efficiency in reporter cell line after tansient transfection
sys.args <- paste("-e",shQuote("rmarkdown::render('mRNA-spliceCount.R')"), "$PWD config_mRNA_Syn9-10_transient.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_mRNA_Syn9-10_transient.txt")
system2("mv", args="mRNA-spliceCount.pdf output/mRNA-spliceCount_mRNA_Syn9-10_transient.pdf")

# Evaluation of RNA-editing efficiency in reporter cell line after stable transduction
sys.args <- paste("-e",shQuote("rmarkdown::render('mRNA-spliceCount.R')"), "$PWD config_mRNA_Syn9-10_stable.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_mRNA_Syn9-10_stable.txt")
system2("mv", args="mRNA-spliceCount.pdf output/mRNA-spliceCount_mRNA_Syn9-10_stable.pdf")

# Evaluation of RNA-editing efficiency in vivo after lentiviral transduction
sys.args <- paste("-e",shQuote("rmarkdown::render('mRNA-spliceCount.R')"), "$PWD config_mRNA_Syn9-10_inVivo.txt", sep=" ")
system2("/usr/local/bin/Rscript", args = sys.args, stdout = "logs/log_mRNA_Syn9-10_inVivo.txt")
system2("mv", args="mRNA-spliceCount.pdf output/mRNA-spliceCount_mRNA_Syn9-10_inVivo.pdf")

#Analysis of Tau Platerunner data
system2("xvfb-run", args = "-a ImageJ-linux64 -macro /home/rstudio/macros/PlaterunnerHeadless.ijm", stdout = "logs/log_TauPlaterunner.txt")
