import subprocess
import argparse
import os

# def get_args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    parser.add_argument("-p", "--profile", 
                        help="path to a snakemake profile, such as a profile to submit jobs to SLURM",
                        default=None)
    parser.add_argument("-c", "--cores", help="Number of cores for snakemake to use",
                        type=int, default=None)
    return parser.parse_args()


def get_snakefile_path():
    path_2_script = os.path.dirname(__file__)
    snakepath = os.path.join(path_2_script, "../snakefile")
    # remove ".."
    return os.path.normpath(snakepath)


def construct_snakemake_command(snakepath, args):
    command = ["snakemake", 
                "-s", snakepath,
                "--configfile", args.config, 
                "--latency-wait", "10", 
                "--rerun-incomplete", 
                "--use-conda"]
    
    if args.profile is not None:
        command.append("--profile")
        command.append(args.profile)

    if args.cores is not None:
        command.append("--cores")
        command.append(str(args.cores))

    print(command)
    return command


def main():
    snakepath = get_snakefile_path()
    args = get_args()

    command = construct_snakemake_command(snakepath, args)
    subprocess.run(command)

if __name__=="__main__":
    main()