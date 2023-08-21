import subprocess
import os
import random
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import chain

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


def sample_from_genome(genome, nreads, microbial):
    # with replacement
    reads = random.choices(genome, k=nreads)
    seqs = [read.seq for read in reads]

    if microbial == True:
        proteins = [get_protein_name_from_fasta_description(read.description) 
                    for read in reads]
    else:
        proteins = [read.description for read in reads]
    return dict(zip(proteins, seqs))


def sample_reads_from_seqs(seqs, read_length):
    for protein, seq in seqs.items():
        # only get a shorter read if the seq is long enough. 
        # Otherwise just leave the seq as is
        if len(seq) > read_length:
            start = random.randint(0,len(seq)-read_length)
            new_seq = seq[start:start+read_length]
            seqs[protein] = new_seq
    return seqs

def get_sampled_reads_from_all_genomes(genomes_read_num_dict, genomes_paths):
    all_reads = [None]*len(genomes_paths)
    for i, genome_name in enumerate(genomes_read_num_dict.keys()):
        genome_path = genomes_paths[i]

        if "human" in genome_name:
            genome = read_human_pangenome(genome_path)
            sampled_seqs = sample_from_genome(genome, genomes_read_num_dict[genome_name], microbial=False)
        else:
            genome = read_microbial_genome(genome_path)
            sampled_seqs = sample_from_genome(genome, genomes_read_num_dict[genome_name], microbial=True)
        
        
        sampled_reads = sample_reads_from_seqs(sampled_seqs, 150)
        all_reads[i] = sampled_reads
    
    return all_reads
        

def simulate_qual_score(mean_phred, var_phred, read_len, min_phred=10):
    # simulate from 40 "trials"
    n = 40
    # setup for simulation from neg binom model
    mean_nb = n - mean_phred
    p = mean_nb/(var_phred**2) # calculate p and r (successes) for n binom
    r = mean_nb**2/(var_phred**2 - mean_nb) # this comes from wikipedia

    # instantiate rng
    rng = np.random.default_rng(seed=42)

    to_subtract = rng.negative_binomial(n=r, p=p, size=read_len)
    # subtract from max phred
    phreds = 40 - to_subtract
    # don't let it be lower than the min phred
    # probs not the best way to do this but good enough?
    phreds = [max(phred, min_phred) for phred in phreds]

    return phreds



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
    sampled_reads_bio = [[SeqRecord(Seq(seq), seq_id, '', '') 
                          for seq_id, seq in genome_sampled_reads.items()] 
                         for genome_sampled_reads in sampled_reads]
    
    sampled_reads_flattened = list(chain.from_iterable(sampled_reads_bio))

    for read in sampled_reads_flattened:
        read.letter_annotations["phred_quality"] = simulate_qual_score(35, 5, len(read.seq))

    print(sampled_reads_flattened)

    SeqIO.write(sampled_reads_flattened, "tests/mock_community/mock_community.fastq", "fastq")


if __name__=="__main__":
    main()