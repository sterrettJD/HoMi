import os
import argparse
import subprocess

def get_args():
    parser = argparse.ArgumentParser(description="Builds a big Bowtie2 reference for human RNA decontamination.")
    parser.add_argument("-o", "--output_dir", type=str, default="host_reference",
                        help="Output directory to store the reference files (default: host_reference)")
    args = parser.parse_args()
    return args


def download_and_decompress(output_dir, name, url):
    """
    Download and decompress a reference file if not already present.
    """
    local_fa = os.path.join(output_dir, f"{name}.fa")

    if url.endswith(".gz"):
        download_path = os.path.join(output_dir, f"{name}.fa.gz")
    else:
        download_path = os.path.join(output_dir, f"{name}.fa")

    if not os.path.exists(local_fa):
        # Skip if already decompressed
        if not os.path.exists(download_path):
            print(f"Downloading {name}...")
            subprocess.run(["wget", "-O", download_path, url], check=True)

        if (url.endswith(".gz")):
            print(f"Decompressing {name}...")
            subprocess.run(["gunzip", "-k", download_path], check=True)


def concatenate_fasta(references, output_filepath, output_dir):
    """
    Concatenate all downloaded reference FASTA files into a single file.
    """
    print("Concatenating reference FASTA files...")
    with open(output_filepath, "w") as outfile:
        for name in references.keys():
            ref_file = os.path.join(output_dir, f"{name}.fa")
            with open(ref_file, "r") as infile:
                outfile.write(infile.read())


def build_bowtie2_index(concatenated_fasta_path, bowtie2_index_prefix):
    """
    Build the Bowtie2 index from the concatenated FASTA file.
    """
    print("Building Bowtie2 index...")
    subprocess.run(["bowtie2-build", " --threads", "4",
                    concatenated_fasta_path, bowtie2_index_prefix], check=True)
    print(f"Bowtie2 index built successfully at: {bowtie2_index_prefix}")


def main():
    args = get_args()
    output_dir = args.output_dir

    # Define URLs for the reference sequences
    REFERENCES = {
        "CHM13v2.0_transcriptome": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_rna.fna.gz",
        "HLA": "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta"
    }

    concatenated_fasta_path = os.path.join(output_dir, "concatenated.fa")
    bowtie2_index_prefix = os.path.join(output_dir, "index")

    os.makedirs(output_dir, exist_ok=True)
    
    for name, url in REFERENCES.items():
        download_and_decompress(output_dir, name, url)


    concatenate_fasta(REFERENCES, concatenated_fasta_path, output_dir)
    build_bowtie2_index(concatenated_fasta_path, bowtie2_index_prefix)


if __name__ == "__main__":
    main()
