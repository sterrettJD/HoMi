import subprocess
import os
import random
from Bio import SeqIO


def fetch_data(accession_id, output_dir, threads=1):
    """
    Pulls data for a given SRA accession using the fasterq-dump.

    Args:
        accession_id (str): The accession number.
        output_dir (str): The directory where the downloaded data will be saved.
        threads (int): The number of threads fasterq-dump should use.

    Returns:
        None
    """
    command = ["fasterq-dump", "-O", output_dir, "--split-files", accession_id, "-e", str(threads)]
    try:
        res = subprocess.run(command, check=True)
        print("Data downloaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output}")
        print("Failed to download data.")


def sample_reads(filepath, n_reads, output_filepath, seed=1234):
    """
    Samples reads from a fastq

    Args:
        filepath (str): The fastq filepath.
        n_reads (int): The number of reads to sample
        output_filepath (str): The file to store the sampled sequences in.

    Returns:
        None


    with open(os.path.join(filepath, FASTQFILEHERE)) as f:
        # should be better to add in some sort of subsetting for sampling without reading all of the file into memory...
        genome = [read for read in SeqIO.FastaIO.FastqIterator(f)]
    random.seed(seed)
    subset = random.choices(genome, k=n_reads)
    
    # NEED TO WRITE TO FASTQ
    """
    
    return None    
    

def main():
    colon_sample_htx="SRP127360"
    output_dir="semisynthetic"
    fetch_data(colon_sample_htx, output_dir)

if __name__ == "__main__":
    main()