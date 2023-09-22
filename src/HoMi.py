import subprocess
import argparse
import os
from re import findall


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
    parser.add_argument("--unlock", action="store_true",
                        help="Pass this to unlock a snakemake directory before running the pipeline.\
                              This is useful if the pipeline failed and you need to rerun it.")
    parser.add_argument("--snakemake_extra",
                        help="Extra parameters to pass snakemake",
                        default=None)
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

    # Maybe just pulling the list of conda envs from that directory could be a better option?
    # But for now, replacing it based on what's in the snakefile seems more direct
    conda_envs_patterns = findall(r"conda_envs\/[a-zA-Z\d\D]+.yaml", file)
    conda_envs_patterns = set(conda_envs_patterns)
    
    print("Conda environments should be named accordingly:")
    for pattern in conda_envs_patterns:
        replacement = pattern
        replacement = replacement.replace("conda_envs/", "")
        replacement = replacement.replace(".yaml", "")

        print(f"{pattern} --> {replacement}")
        file = file.replace(pattern, replacement)


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
    
    if args.snakemake_extra is not None:
        command.extend(args.snakemake_extra.split())

    return command


def main():
    args = get_args()

    snakepath = get_snakefile_path()

    if args.conda_prebuilt == True:
        snakepath = make_prebuilt_conda_snakefile(snakepath)

    if args.unlock == True:
        subprocess.run(["snakemake", 
                        "-s", snakepath,
                        "--configfile", args.config,
                        "--unlock"])
    
    command = construct_snakemake_command(snakepath, args)
    completed_process = subprocess.run(command)

if __name__=="__main__":
    main()
