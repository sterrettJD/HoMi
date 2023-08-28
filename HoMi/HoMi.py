import subprocess
def main():
    subprocess.run(["snakemake", 
                    "--configfile", "tests/example_config.yaml", 
                    "-c1", 
                    "--latency-wait", "10", 
                    "--rerun-incomplete", "--use-conda"])

if __name__=="__main__":
    main()