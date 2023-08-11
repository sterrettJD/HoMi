import subprocess

def main():
    subprocess.run(["bash", "tests/download_genomes_for_mock.sh"])

if __name__=="__main__":
    main()