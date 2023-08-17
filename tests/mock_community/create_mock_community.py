import subprocess
import os

def download_genomes(path="mock_community/data"):
    subprocess.run(["cd", path])
    subprocess.run(["bash", "tests/download_genomes_for_mock.sh"])


def main():
    if os.path.exists("mock_community/data")==F:
        download_genomes(path="mock_community/data")
    # TODO: add a check to make sure it has the right stuff

if __name__=="__main__":
    main()