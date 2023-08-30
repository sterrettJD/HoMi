import subprocess
import argparse
import os

# def get_args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    return parser.parse_args()


def get_snakefile_path():
    path_2_script = os.path.dirname(__file__)
    snakepath = os.path.join(path_2_script, "../snakefile")
    # remove ".."
    return os.path.normpath(snakepath)


def main():
    snakepath = get_snakefile_path()

    args = get_args()

    subprocess.run(["snakemake", 
                    "-s", snakepath,
                    "--configfile", args.config, 
                    "-c1", 
                    "--latency-wait", "10", 
                    "--rerun-incomplete", 
                    "--use-conda"])

if __name__=="__main__":
    main()

# https://stackoverflow.com/questions/74260602/how-to-run-a-python-script-relative-to-its-location