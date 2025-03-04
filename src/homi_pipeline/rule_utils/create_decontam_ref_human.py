import os
import subprocess

def main():

    # Define URLs for the reference sequences (update if needed)
    references = {
        "CHM13v2.0": "https://example.com/path/to/CHM13v2.0.fa.gz",
        "GRCh38": "https://example.com/path/to/GRCh38.fa.gz",
        "HLA": "https://example.com/path/to/HLA.fa.gz",
        "GRCh38_transcriptome": "https://example.com/path/to/GRCh38_transcriptome.fa.gz"
    }

    output_dir = "host_reference"
    concatenated_fasta = os.path.join(output_dir, "combined_reference.fa")
    bowtie2_index_prefix = os.path.join(output_dir, "host_index")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Download and decompress reference files
    for name, url in references.items():
        local_file = os.path.join(output_dir, f"{name}.fa.gz")
        decompressed_file = local_file.rstrip(".gz")

        if not os.path.exists(local_file):  # Avoid re-downloading
            subprocess.run(["wget", "-O", local_file, url], check=True)

        # Decompress if necessary
        if not os.path.exists(decompressed_file):
            subprocess.run(["gunzip", "-k", local_file], check=True)

    # Concatenate all reference files
    with open(concatenated_fasta, "w") as outfile:
        for name in references.keys():
            ref_file = os.path.join(output_dir, f"{name}.fa")
            with open(ref_file, "r") as infile:
                outfile.write(infile.read())

    # Build Bowtie2 index
    subprocess.run(["bowtie2-build", concatenated_fasta, bowtie2_index_prefix], check=True)

    print(f"Bowtie2 index built successfully at: {bowtie2_index_prefix}")

if __name__=="__main__":
    main()