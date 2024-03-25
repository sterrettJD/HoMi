import subprocess

def fetch_data(accession_id, output_dir, threads=1):
    """
    Pulls data for a given SRA accession using the fasterq-dump.

    Args:
        accession_id (str): The accession number.
        output_dir (str): The directory where the downloaded data will be saved.
        threads (int): The number of threads fasterq-dump should use.

    Returns:
        None
    """
    command = ["fasterq-dump", "-O", output_dir, "--split-files", accession_id, "-e", str(threads)]
    try:
        res = subprocess.run(command, check=True)
        print("Data downloaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output}")
        print("Failed to download data.")


def main():
    colon_sample_htx="SRP127360"
    output_dir="semisynthetic"
    fetch_data(colon_sample_htx, output_dir)

if __name__ == "__main__":
    main()