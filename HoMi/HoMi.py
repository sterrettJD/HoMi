import subprocess
import argparse
import os


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    parser.add_argument("-p", "--profile", 
                        help="path to a snakemake profile, such as a profile to submit jobs to SLURM",
                        default=None)
    parser.add_argument("-c", "--cores", help="Number of cores for snakemake to use",
                        type=int, default=None)
    parser.add_argument("--conda_prebuilt", 
                        action="store_true",
                        help="All conda environments have been prebuilt, snakemake should not build them")
    return parser.parse_args()


def get_snakefile_path():
    path_2_script = os.path.dirname(__file__)
    snakepath = os.path.join(path_2_script, "../snakefile")
    # remove ".."
    return os.path.normpath(snakepath)


def make_prebuilt_conda_snakefile_name(old_snakepath):
    return old_snakepath + "_prebuilt"


def make_prebuilt_conda_snakefile(snakepath):
    with open(snakepath, "r") as f:
        file = f.read()

    file = file.replace("conda_envs/", "")

    new_snakepath = make_prebuilt_conda_snakefile_name(snakepath)
    with open(new_snakepath, "w") as f:
        f.write(file)

    return new_snakepath


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

    return command


def main():
    args = get_args()

    snakepath = get_snakefile_path()

    if args.conda_prebuilt is True:
        snakepath = make_prebuilt_conda_snakefile(snakepath)

    command = construct_snakemake_command(snakepath, args)
    subprocess.run(command)

if __name__=="__main__":
    main()