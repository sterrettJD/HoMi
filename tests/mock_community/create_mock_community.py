import subprocess
import os
import random
from Bio import SeqIO

def download_genomes(fpath=os.path.join("tests", "mock_community" ,"data")):
    print(os.getcwd())
    os.chdir(fpath)
    # make sure we're in the right place
    print(f"Downloading genomes in: {os.getcwd()}")

    subprocess.run(["bash", os.path.join("..", "download_genomes_for_mock.sh")])
    # maybe return a list of the genomes?


def read_microbial_genome(filepath):

    with open(os.path.join(filepath, "cds_from_genomic.fna")) as f:
        genome = [read for read in SeqIO.FastaIO.FastaIterator(f)]
    return genome


def get_protein_name_from_fasta_description(description):
    protein_with_extra = description.split("protein=")[-1]
    protein = protein_with_extra.split("]")[0]
    return protein


def sample_from_genome(genome, nreads):
    # with replacement
    reads = random.choices(genome, k=nreads)
    seqs = [read.seq for read in reads]

    proteins = [get_protein_name_from_fasta_description(read.description) 
                for read in reads]
    return dict(zip(proteins, seqs))

    


def main():
    if os.path.exists("tests/mock_community/data") == False:
        os.mkdir("tests/mock_community/data/")
        print("Creating directory tests/mock_community/data")

    genomes = ["e_coli", "c_beijerinckii", "f_prausnitzii", "human_pangenome"]
    genomes_paths = [os.path.join('tests', 'mock_community', 'data', g, 'genome')
                     for g in genomes]
    genomes_exist = [os.path.exists(genome_path) 
                     for genome_path in genomes_paths]

    if sum(genomes_exist) < len(genomes):
       print("genomes needed")
       # download_genomes()
    else:
        print("All genomes have already been downloaded.")
    # TODO: add a check to make sure it has the right stuff?

    # sample some reads with replacement
    random.seed(42)
    genome = read_microbial_genome(genomes_paths[0])
    sampled_reads = sample_from_genome(genome, 2)
    print(sampled_reads)

if __name__=="__main__":
    main()