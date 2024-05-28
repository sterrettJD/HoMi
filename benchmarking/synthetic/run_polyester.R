#!/usr/bin/env Rscript

# example usage
# Rscript run_polyester.R -t data/host_transcriptome.fna.gz \
# -tu https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
# -g data/host_transcriptome.gff.gz \
# -gu https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
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
    return (unlist(vals))
}





#### MAIN #####
option_list <- list( 
    make_option(c("-t", "--transcriptome_filepath"),
                help="Reference transcriptome filepath"),
    make_option(c("--transcriptome_url"), 
                help="URL to download the reference transcriptome from if it isn't present."),
    make_option(c("-g", "--gtf_filepath"), 
                help="Reference gtf filepath"),
    make_option(c("--gtf_url"), 
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
transcriptome_filepath <- opt$transcriptome_filepath
transcriptome_url <- opt$transcriptome_url
gtf_filepath <- opt$gtf_filepath
gtf_url <- opt$gtf_url
sample_data <- opt$sample_data
host_id <- opt$host_id

get_reference(transcriptome_filepath, transcriptome_url)
get_reference(gtf_filepath, gtf_url)


depths <- get_sample_depths(sample_data, host_id)
nonzero_depths <- depths[depths > 0]
smallest_depth <- min(nonzero_depths)
fold_changes <- nonzero_depths/smallest_depth
groups <- unique(fold_changes)
size_per_group <- as.numeric(table(as.factor(nonzero_depths)))

numtx <- length(readDNAStringSet(transcriptome_filepath))

fold_change_values <- sapply(groups, function(group) rep(group, numtx))
fold_change_matrix <- matrix(fold_change_values, nrow=numtx, byrow=FALSE)


simulate_experiment(fasta=transcriptome_filepath,
                    reads_per_transcript=rnbinom(n=numtx, size=20, prob=0.5),
                    fold_changes=fold_change_matrix,
                    outdir='simulated_reads',
                    num_reps=size_per_group)