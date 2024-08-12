import os
import subprocess
import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_data", help="A CSV file with the columns 'genome', 'GCF_id',"
                                            "and a column for each sample that should be created. "
                                            "Rows denote genomes to sample reads from, and the values "
                                            "in each sample's column denotes how many reads should come "
                                            "from each genome. At least one genome should be 'human', "
                                            "and the GCF id for it does not matter.")
    parser.add_argument("--work_dir", help="The working directory in which data should be downloaded and "
                                            "The <sample_data> path should NOT be relative to this path "
                                            "(it should be relative to where the) script is run from.")
    return parser.parse_args()

def parse_sample_data(filepath):
    df = pd.read_csv(filepath)
    return df

def download_microbial_genome(genome_name, accession_id, fpath):
    start_dir = os.getcwd()
    print(f"Moving out of {start_dir}")
    os.chdir(fpath)
    # make sure we're in the right place
    print(f"Downloading genomes in: {os.getcwd()}")

    # setup
    os.makedirs(genome_name, exist_ok=True)
    os.chdir(genome_name)

    # Download and unzip
    subprocess.run(["curl", "-OJX", "GET", f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession_id}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename={accession_id}.zip", "-H", "Accept: application/zip"])
    subprocess.run(["unzip", f"{accession_id}.zip"])
    
    # cleanup
    subprocess.run(["rm", f"{accession_id}.zip"])
    subprocess.run(["mv", f"ncbi_dataset/data/{accession_id}", "./genome"])
    os.chdir(start_dir)


def download_human_pangenome(fpath):
    start_dir = os.getcwd()
    print(f"Moving out of {start_dir}")
    os.chdir(fpath)
    # make sure we're in the right place
    print(f"Downloading human pangenome in: {os.getcwd()}")

    # setup
    os.makedirs("human", exist_ok=True)
    os.chdir("human")
    
    subprocess.run(["wget", "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/CHM13v11Y.fa.gz"])
    subprocess.run(["gunzip", "CHM13v11Y.fa.gz"])
    
    os.makedirs("genome", exist_ok=True)
    subprocess.run(["mv", "CHM13v11Y.fa", "genome/"])
    os.chdir(start_dir)


def main():
    # get arguments
    args = get_args()
    sample_data = parse_sample_data(args.sample_data)

    # prep working directory and output
    work_dir = args.work_dir
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)

    if os.path.exists("data") == False:
        os.mkdir("data/")
        print(f"Creating directory {work_dir}/data")


    # list of genomes to include
    genomes = sample_data["genome"].to_list()
    accession_ids = sample_data["GCF_id"].to_list()

    # Get non-host genomes and accession IDs
    nh_genomes = [x for x in genomes if "human" not in x]
    nh_accession_ids = [accession_ids[i] for i, x in enumerate(genomes) if "human" not in x]
    
    genomes_paths = [os.path.join("data", g, "genome")
                     for g in genomes]
    # check if these genomes exist
    genomes_exist = [os.path.exists(genome_path) 
                     for genome_path in genomes_paths]
    
    if sum(genomes_exist) < len(genomes):
        print("At least one genome is missing. "
             "Redownloading all genomes, as others might be incomplete or missing.")
        for i, genome in enumerate(nh_genomes):
           download_microbial_genome(genome, nh_accession_ids[i], fpath="data")
        download_human_pangenome(fpath="data")
    else:
        print("All genomes have already been downloaded.")


if __name__=="__main__":
    main()
