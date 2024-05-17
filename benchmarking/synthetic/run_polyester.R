#!/usr/bin/env Rscript

# example usage
# Rscript run_polyester.R -t data/host_transcriptome.fna.gz \
# -u https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
# -s sample_data.csv \
# -n human


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
    message(paste("File downloaded from", url, "to", file_path))
    } else {
    message("File already exists: ", file_path)
    }
}


get_sample_depths <- function(sample_data, host_id){
    df <- read.csv(sample_data)
    rownames(df) <- df$genome
    df$genome <- NULL
    
    index <- !grepl(pattern="GCF_id", x=colnames(df))
    vals <- df[host_id, index]
    return (vals)
}


#### MAIN #####
option_list <- list( 
    make_option(c("-t", "--transcriptome_filepath"),
                help="Reference transcriptome filepath"),
    make_option(c("-u", "--transcriptome_url"), 
                help="URL to download the reference transcriptome from if it isn't present."),
    make_option(c("-s", "--sample_data"), 
                help="A CSV file with the number of reads per sample belonging to the host. Samples are columns, and rows are organisms. A column labeled `genome` should denote the organisms, and a column `GCF_id` will be ignored."),
    make_option(c("-n", "--host_id"), 
                help="A CSV file with the number of reads per sample belonging to the host.")

)

# set this so the download doesn't time out
options(timeout=max(300, getOption("timeout")))

opt <- parse_args(OptionParser(option_list=option_list))

# Extract option values
file_path <- opt$transcriptome_filepath
url <- opt$transcriptome_url
sample_data <- opt$sample_data
host_id <- opt$host_id

get_reference(file_path, url)
depths <- get_sample_depths(sample_data, host_id)
print(unlist(depths))