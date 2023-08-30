import subprocess
import os

# def get_args

def get_snakefile_path():
    path_2_script = os.path.dirname(__file__)
    snakepath = os.path.join(path_2_script, "../snakefile")
    # remove ".."
    return os.path.normpath(snakepath)


def main():
    snakepath = get_snakefile_path()
    subprocess.run(["snakemake", 
                    "-s", snakepath,
                    "--configfile", "tests/example_config.yaml", 
                    "-c1", 
                    "--latency-wait", "10", 
                    "--rerun-incomplete", 
                    "--use-conda"])

if __name__=="__main__":
    main()

# https://stackoverflow.com/questions/74260602/how-to-run-a-python-script-relative-to-its-location