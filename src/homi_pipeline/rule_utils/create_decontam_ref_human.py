import os
import argparse
import subprocess

# Define URLs for the reference sequences
REFERENCES = {
    "CHM13v2.0_transcriptome": "https://example.com/path/to/CHM13v2.0.fa.gz",
    "GRCh38": "https://example.com/path/to/GRCh38.fa.gz",
    "HLA": "https://example.com/path/to/HLA.fa.gz",
    "GRCh38_transcriptome": "https://example.com/path/to/GRCh38_transcriptome.fa.gz"
}

OUTPUT_DIR = "host_reference"
CONCATENATED_FASTA = os.path.join(OUTPUT_DIR, "combined_reference.fa")
BOWTIE2_INDEX_PREFIX = os.path.join(OUTPUT_DIR, "host_index")


def get_args():
    parser = argparse.ArgumentParser(description="Builds a big Bowtie2 reference for human RNA decontamination.")
    parser.add_argument("-o", "--output_dir", type=str, default="host_reference",
                        help="Output directory to store the reference files (default: host_reference)")
    args = parser.parse_args()
    return args


def download_and_decompress(name, url):
    """
    Download and decompress a reference file if not already present.
    """
    local_gz = os.path.join(OUTPUT_DIR, f"{name}.fa.gz")
    local_fa = os.path.join(OUTPUT_DIR, f"{name}.fa")

    if not os.path.exists(local_fa):  # Skip if already decompressed
        if not os.path.exists(local_gz):
            print(f"Downloading {name}...")
            subprocess.run(["wget", "-O", local_gz, url], check=True)

        print(f"Decompressing {name}...")
        subprocess.run(["gunzip", "-k", local_gz], check=True)


def concatenate_fasta():
    """
    Concatenate all downloaded reference FASTA files into a single file.
    """
    print("Concatenating reference FASTA files...")
    with open(CONCATENATED_FASTA, "w") as outfile:
        for name in REFERENCES.keys():
            ref_file = os.path.join(OUTPUT_DIR, f"{name}.fa")
            with open(ref_file, "r") as infile:
                outfile.write(infile.read())


def build_bowtie2_index():
    """
    Build the Bowtie2 index from the concatenated FASTA file.
    """
    print("Building Bowtie2 index...")
    subprocess.run(["bowtie2-build", CONCATENATED_FASTA, BOWTIE2_INDEX_PREFIX], check=True)
    print(f"Bowtie2 index built successfully at: {BOWTIE2_INDEX_PREFIX}")


def main():
    args = get_args()
    output_dir = args.output_dir

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    for name, url in REFERENCES.items():
        download_and_decompress(name, url)

    concatenate_fasta()
    build_bowtie2_index()


if __name__ == "__main__":
    main()
