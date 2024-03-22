import subprocess

def fetch_data(accession_id, output_dir):
    """
    Pulls data for a given BioProject accession using the ncbi-datasets tool.

    Args:
        accession_id (str): The accession number.
        output_dir (str): The directory where the downloaded data will be saved.

    Returns:
        None
    """
    command = ["fastq-dump", "-O", output_dir, "--split-files", accession_id]
    try:
        res = subprocess.run(command, check=True)
        print("Data downloaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output}")
        print("Failed to download data.")


def main():
    colon_sample_htx="SRR28258513"
    output_dir="semisynthetic"
    fetch_data(colon_sample_htx, output_dir)

if __name__ == "__main__":
    main()