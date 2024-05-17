#!/usr/bin/env Rscript

if(!require("optparse")){
    install.packages("optparse")
}

if (!requireNamespace("BiocManager", quietly=TRUE)){
    install.packages("BiocManager")
}

if (!require("polyester")){
    BiocManager::install("polyester")
}

if (!require("Biostrings")){
    BiocManager::install("Biostrings")
}

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("polyester"))
suppressPackageStartupMessages(library("Biostrings"))


get_reference <- function(file_path, url){
    # Check if the file exists
    if (!file.exists(file_path)) {
    # If file doesn't exist, download from URL
    download.file(url=url, destfile=file_path, mode = "wb")
    message(paste("File downloaded from", download_url, "to", file_path))
    } else {
    message("File already exists:", file_path)
    }
}




#### MAIN #####
option_list <- list( 
    make_option(c("-t", "--transcriptome_filepath"),
                help="Reference transcriptome filepath"),
    make_option(c("-u", "--transcriptome_url"), 
                help="URL to download the reference transcriptome from if it isn't present."),
    make_option(c("-s", "--sample_data"), 
                help="A CSV file with the number of reads per sample belonging to the host.")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Extract option values
file_path <- opt$transcriptome_filepath
url <- opt$transcriptome_url

get_reference(file_path, url)