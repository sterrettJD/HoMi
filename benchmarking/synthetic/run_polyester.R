#!/usr/bin/env Rscript

# example usage
# Rscript run_polyester.R -t data/host_transcriptome.fna.gz \
# --transcriptome_url https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
# -g data/host_transcriptome.gff.gz \
# --gtf_url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
# -s sample_data.csv \
# -n human \
# -o synthetic_transcriptomes


if(!require("optparse")){
    install.packages("optparse", repos="http://cran.us.r-project.org")
}

if (!requireNamespace("BiocManager", quietly=TRUE)){
    install.packages("BiocManager", repos="http://cran.us.r-project.org")
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
    df <- read.csv(sample_data, check.names=FALSE)
    rownames(df) <- df$genome
    df$genome <- NULL

    index <- !grepl(pattern="GCF_id", x=colnames(df))
    vals <- df[host_id, index]
    return (unlist(vals))
}


check_reads_per_transcript <- function(reads, numtx){
  if (reads < numtx){
    warning("Number of reads per transcript is less than 1 in the baseline group. \
            Please consider scaling up.")
  }
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
                help="A CSV file with the number of reads per sample belonging to the host."),
    make_option(c("-o", "--output_dir"), 
                help="The output directory for reads.")

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
output_dir <- opt$output_dir

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


read_per_transcript <- smallest_depth/numtx
check_reads_per_transcript(smallest_depth, numtx)
print(paste("Reads per transcript (baseline group):", read_per_transcript))
print(paste("Total reads (baseline group):", sum(rep(read_per_transcript, numtx))))
print("Fold change matrix (head):")
print(head(fold_change_matrix))

simulate_experiment(fasta=transcriptome_filepath,
                    reads_per_transcript=rep(read_per_transcript, numtx),
                    fold_changes=fold_change_matrix,
                    outdir=output_dir,
                    num_reps=size_per_group,
                    readlen=150)

# rename samples according to the provided IDs
start_wd <- getwd()
setwd(output_dir)
i <- 1
for(sample in names(nonzero_depths)){
    if (i < 10){
        s1 <- paste0("sample_0", i, "_1.fasta")
        s2 <- paste0("sample_0", i, "_2.fasta")
    } else {
        s1 <- paste0("sample_", i, "_1.fasta")
        s2 <- paste0("sample_", i, "_2.fasta")
    }
    print(paste(s1, "->", sample))

    file.rename(s1, paste0(sample, "_R1.fasta"))
    file.rename(s2, paste0(sample, "_R2.fasta"))

    i <- i + 1
}

# create the 0 reads files
zero_depths <- depths[depths == 0]
for(sample in names(zero_depths)){
    file.create(paste0(sample, "_R1.fasta"))
    file.create(paste0(sample, "_R2.fasta"))
}

setwd(start_wd)
