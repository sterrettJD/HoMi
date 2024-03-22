import subprocess

def fetch_data(bioproj, output_dir):
    """
    Pulls data for a given BioProject accession using the ncbi-datasets tool.

    Args:
        bioproj (str): The BioProject accession number.
        output_dir (str): The directory where the downloaded data will be saved.

    Returns:
        None
    """
    command = ["datasets", "download", "genome", "accession", "--id", bioproj, "--out", output_dir]
    try:
        res = subprocess.run(command, check=True)
        print("Data downloaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output}")
        print("Failed to download data.")


def main():
    HMP2_htx="PRJNA438663"
    output_dir="semisynthetic"
    fetch_data(HMP2_htx, output_dir)

if __name__ == "__main__":
    main()