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


def read_human_pangenome(filepath):
    with open(os.path.join(filepath, "CHM13v11Y.fa")) as f:
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


def get_sampled_reads_from_all_genomes(genomes_read_num_dict, genomes_paths):
    all_reads = [None]*len(genomes_paths)
    for genome_name, i in enumerate(genomes_read_num_dict):
        genome_path = genomes_paths[i]

        if genome_name.contains("human"):
            genome = read_human_pangenome(genome_path)
        else:
            genome = read_microbial_genome(genome_path)
        
        sampled_reads = sample_from_genome(genome, genomes_read_num_dict[genome_name])
        all_reads[i] = sampled_reads
    
    return all_reads
        


def main():
    # check for mock community data already existing
    if os.path.exists("tests/mock_community/data") == False:
        os.mkdir("tests/mock_community/data/")
        print("Creating directory tests/mock_community/data")

    # list of genomes to include
    genomes = ["e_coli", "c_beijerinckii", "f_prausnitzii", "human_pangenome"]
    genomes_paths = [os.path.join('tests', 'mock_community', 'data', g, 'genome')
                     for g in genomes]
    
    # check if these genomes exist
    genomes_exist = [os.path.exists(genome_path) 
                     for genome_path in genomes_paths]
    if sum(genomes_exist) < len(genomes):
       print("At least one genome is missing. "
             "Redownloading all genomes, as others might be incomplete or missing.")
       download_genomes()
    else:
        print("All genomes have already been downloaded.")
    # TODO: add a check to make sure it has the right stuff?

    # sample some reads with replacement
    random.seed(42)
    genomes_read_num_dict = {
        "e_coli": 50, 
        "c_beijerinckii": 50, 
        "f_prausnitzii": 50, 
        "human_pangenome": 150
    }
    sampled_reads = get_sampled_reads_from_all_genomes(genomes_read_num_dict, genomes_paths)


    

if __name__=="__main__":
    main()