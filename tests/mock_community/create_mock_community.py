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


def subsample_reads():
    random.seed(42)
    # do things


def main():
    if os.path.exists("tests/mock_community/data") == False:
        os.mkdir("tests/mock_community/data/")
        print("Creating directory tests/mock_community/data")

    genomes = ["e_coli", "c_beijerinckii", "f_prausnitzii", "human_pangenome"]
    genomes_exist = [os.path.exists(os.path.join('tests', 'mock_community', 'data', g)) 
                     for g in genomes]

    if sum(genomes_exist) < len(genomes):
        download_genomes()
    else:
        print("All genomes have already been downloaded.")
    # TODO: add a check to make sure it has the right stuff?


if __name__=="__main__":
    main()