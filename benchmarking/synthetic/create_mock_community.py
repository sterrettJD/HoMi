import subprocess
import os
import random
import numpy as np
import pandas as pd
import gzip
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import chain


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_data", help="A CSV file with the columns 'genome', 'GCF_id',"
                                            "and a column for each sample that should be created. "
                                            "Rows denote genomes to sample reads from, and the values "
                                            "in each sample's column denotes how many reads should come "
                                            "from each genome. At least one genome should be 'human', "
                                            "and the GCF id for it does not matter.")
    parser.add_argument("sample_id", help="The sample ID that should be created. "
                                          "This should be a column in the sample_data file.")
    parser.add_argument("--work_dir", help="The working directory relative to which the data should be downloaded and "
                                            "fastq files should be created. The <sample_data> path should "
                                            "NOT be relative to this path (it should be relative to where the"
                                            "script is run from).")
    parser.add_argument("--leave_unzipped", help="Pass this flag if unzipped fastq files should be left. "
                                            "Otherwise, they will be deleted",
                        action="store_true")
    parser.add_argument("--output_dir", help="Directory (relative to work_dir) where output fastqs should be written.")
    return parser.parse_args()


def parse_sample_data(filepath):
    df = pd.read_csv(filepath)
    return df


def check_sample_id_exists(sample_id, samples):
    if sample_id not in samples:
        raise ValueError(f"{sample_id} does not exist as a column in the sample data.")


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

    # Currently not using protein name but can maybe use this later
    if microbial == True:
        proteins = [get_protein_name_from_fasta_description(read.description) 
                    for read in reads]
    else:
        proteins = [read.description for read in reads]
    return dict(zip(np.arange(len(seqs)), seqs))


def sample_reads_from_seqs(seqs, read_length):
    for i, seq in seqs.items():
        # only get a shorter read if the seq is long enough. 
        # Otherwise just leave the seq as is
        if len(seq) > read_length:
            start = random.randint(0,len(seq)-read_length)
            new_seq = seq[start:start+read_length]
            seqs[i] = new_seq
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
        

def simulate_qual_score(mean_phred, var_phred, read_len, min_phred=10, method="gaussian"):
    # instantiate rng
    rng = np.random.default_rng(seed=random.randint(0,9999))

    if method=="gaussian":
        phreds = rng.normal(mean_phred, var_phred, size=read_len)
        phreds = phreds.astype(int)

    elif method=="negative_binomial":
        # simulate from 40 "trials"
        n = 40
        # setup for simulation from neg binom model
        mean_nb = n - mean_phred
        p = mean_nb/(var_phred**2) # calculate p and r (successes) for n binom
        r = mean_nb**2/(var_phred**2 - mean_nb) # this comes from wikipedia

        to_subtract = rng.negative_binomial(n=r, p=p, size=read_len)

        # subtract from max phred
        phreds = 40 - to_subtract
    
    else:
        raise ValueError(f"argument method {method} is not available. Please use 'gaussian' or 'negative_binomial'")

    # don't let it be lower than the min phred
    # probs not the best way to do this but good enough?
    phreds = [max(phred, min_phred) for phred in phreds]

    return phreds

def get_error_probs(phred_score):
    # just using the PHRED score definition here
    error = 1/(10**(phred_score/10))
    no_error = 1-error
    return[no_error, error]


def mutate_seq_from_qual(seqrecord):
    letters = {"A", "C", "T", "G"}
    # decide whether or not to mutate each base pair
    to_mut = [random.choices([0,1], weights=get_error_probs(q))[0] 
              for q in seqrecord.letter_annotations["phred_quality"]]
    # gotta convert to an iterable because Seq doesn't allow assignment
    sequence = list(seqrecord.seq)
    for i, base in enumerate(sequence):
        if to_mut[i]==1:
            # naive substitution
            sequence[i] = random.choice(list(letters.difference({base})))
    # convert to string then Seq object
    new_seqrecord = seqrecord
    new_seqrecord.seq = Seq("".join(sequence))

    return new_seqrecord


def create_qual_scores_and_mutate(list_of_seqrecords, 
                                  mean_phred, var_phred, min_phred=10):
    # generate quality scores (somewhat naively)
    for read in list_of_seqrecords:
        read.letter_annotations["phred_quality"] = simulate_qual_score(mean_phred, var_phred, len(read.seq), min_phred)
    
    # mutate the reads based on quality scores
    sampled_reads_flat_mut = [mutate_seq_from_qual(seqrecord) 
                               for seqrecord in list_of_seqrecords]
    return sampled_reads_flat_mut


def compress_fastq(fp, remove_unzipped=True):
    with open(fp, 'rb') as f_in, gzip.open(f"{fp}.gz", 'wb') as f_out:
        f_out.writelines(f_in)

    if remove_unzipped:
        os.remove(fp)


def main():
    # get arguments
    args = get_args()
    sample_data = parse_sample_data(args.sample_data)

    # prep working directory and output
    work_dir = args.work_dir
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)

    if args.output_dir is None:
        output_dir = "."
    else:
        output_dir = args.output_dir
        os.makedirs(output_dir, exist_ok=True)

    # check for mock community data already existing
    if os.path.exists("data") == False:
        os.mkdir("data/")
        print(f"Creating directory {work_dir}/data")

    # list of genomes to include
    genomes = sample_data["genome"].to_list()
    genomes_paths = [os.path.join("data", g, "genome")
                     for g in genomes]
    # check if these genomes exist
    genomes_exist = [os.path.exists(genome_path) 
                     for genome_path in genomes_paths]
    
    if sum(genomes_exist) < len(genomes):
        raise ValueError("At least one genome is missing. "
             "Please redownload all genomes, as others might be incomplete or missing.")
    else:
        print("All genomes have already been downloaded.")

    # Get the sample ID
    sample_columns = [col for col in sample_data.columns if (any(s in col for s in ["genome", "GCF_id"])==False)]
    sample_id = args.sample_id
    check_sample_id_exists(sample_id, sample_columns)
    # Seed is set as a hash of the sample ID, which should be reproducible but is not 
    # identical for each sample
    random.seed(hash(sample_id))

    # create a dictionary of the number of reads from each genome
    genomes_read_num_dict = dict(zip(sample_data["genome"], sample_data[sample_id].astype(int)))
    print(genomes_read_num_dict)
    # Sample the reads from all genomes
    sampled_reads = get_sampled_reads_from_all_genomes(genomes_read_num_dict, genomes_paths)
    sampled_reads_bio = [[SeqRecord(Seq(seq), '', '', '') 
                        for seq_id, seq in genome_sampled_reads.items()] 
                        for genome_sampled_reads in sampled_reads]
    
    # flatten the list of lists (reads per genome) to just a list (reads)
    sampled_reads_flat = list(chain.from_iterable(sampled_reads_bio))

    # create the reverse reads
    sampled_reads_flat_rev = [SeqRecord(Seq(''), '', '', '')] * len(sampled_reads_flat)
    for i, sequence in enumerate(sampled_reads_flat):
        sampled_reads_flat_rev[i].seq = sequence.seq[::-1] 
    
    # Create quality scores and mutate the reads
    sampled_reads_flat_mut = create_qual_scores_and_mutate(sampled_reads_flat, 
                                        mean_phred=35, var_phred=5, min_phred=10)
    sampled_reads_flat_rev_mut = create_qual_scores_and_mutate(sampled_reads_flat_rev, 
                                        mean_phred=35, var_phred=3, min_phred=10)
    
    # Write the reads
    SeqIO.write(sampled_reads_flat_mut, os.path.join(output_dir, f"{sample_id}_R1.fastq"), "fastq")
    SeqIO.write(sampled_reads_flat_rev_mut, os.path.join(output_dir, f"{sample_id}_R2.fastq"), "fastq")

    # Zip the reads
    remove_unzipped = args.leave_unzipped == False
    compress_fastq(os.path.join(output_dir, f"{sample_id}_R1.fastq"), remove_unzipped=remove_unzipped)
    compress_fastq(os.path.join(output_dir, f"{sample_id}_R2.fastq"), remove_unzipped=remove_unzipped)


if __name__=="__main__":
    main()