#!/usr/bin/env python3

import argparse
import os
from cookiecutter.main import cookiecutter


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster_type", 
                        choices=["slurm-smt", "slurm-nosmt"],
                        help="The type of cluster for which the profile should be created. \
                              Currently supported options are slurm-smt (slurm with simultaneous multithreading) \
                              and slurm-nosmt (slurm-managed cluster without simultaneous multithreading).",
                        required=True)
    parser.add_argument("--output_dir", 
                        help="Where should the profile be stored?",
                        default=".")
    return parser.parse_args()


def get_cookiecutter_slurm_profile(output_dir):
    cookiecutter("gh:Snakemake-Profiles/slurm", output_dir=output_dir, no_input=True)


def check_profile_named_slurm(output_dir):
    alleged_profile_dir= os.path.join(output_dir, "slurm")
    if not os.path.isdir(alleged_profile_dir):
        raise ValueError("""
                         The profile's name has not been set to slurm.
                         Please make sure the profile is named slurm
                         when following the prompts for profile setup.
                         """)

    files_should_be_present = ["config.yaml",
                               "CookieCutter.py",
                               "slurm_utils.py",
                               "slurm-jobscript.sh",
                               "slurm-sidecar.py",
                               "slurm-status.py",
                               "slurm-submit.py"]
    contents_of_profile_dir = os.listdir(alleged_profile_dir)
    
    files_there = [file in contents_of_profile_dir 
                   for file in files_should_be_present]
    for i, file in enumerate(files_there):
        if not file:
            raise ValueError(f"""
                             The slurm profile doesn't seem to have all necessary files.
                             {files_should_be_present[i]} is missing.
                             Something may have gone wrong, please check this directory.
                             """)


def write_new_config(filepath, new_contents):
    with open(filepath, 'w') as file:
        file.write(new_contents)


def set_use_conda_true(config_contents):
    return config_contents.replace("use-conda: \"False\"", 
                                   "use-conda: \"True\"")
    

def set_print_shell_true(config_contents):
    return config_contents.replace("printshellcmds: \"False\"", 
                                   "printshellcmds: \"True\"")


def update_HoMi_reqs_in_config(output_dir):
    config_path = os.path.join(output_dir, "slurm", "config.yaml")
    with open(config_path, 'r') as file:
        content = file.read()

    with_conda = set_use_conda_true(content)
    with_conda_and_shell = set_print_shell_true(with_conda)

    write_new_config(config_path, with_conda_and_shell)


def check_use_conda_slurm(output_dir):
    config_path = os.path.join(output_dir, "slurm", "config.yaml")

    with open(config_path, 'r') as file:
        content = file.read()
        if "use-conda: \"True\"" not in content:
            raise ValueError("""
                             Conda is not enabled with your Slurm cluster profile.
                             Please set use_conda to True.
                             """)
        

def write_new_utils(filepath, new_contents):
    with open(filepath, 'w') as file:
        file.write(new_contents)


def convert_slurm_profile_tasks(output_dir):
    slurm_utils_path = os.path.join(output_dir, "slurm", "slurm_utils.py")
    with open(slurm_utils_path, 'r') as file:
        content = file.read()
        # Make sure there is something to replace
        if "options[\"cpus-per-task\"] = job_properties[\"threads\"]" not in content:
            raise ValueError(f"""
                             {slurm_utils_path} does not contain the substring to be replaced to integrate with no SMT. 
                             Please re-evaluate what's going on here...
                             """)
        else:
            # the cpus-per-task param doesn't work when threads >1 and the core isn't hyperthreaded
            # this replaces the line in slurm_utils.py to use ntasks instead of cpus-per-task
            new = content.replace("options[\"cpus-per-task\"] = job_properties[\"threads\"]", 
                                  "options[\"ntasks\"] = job_properties[\"threads\"]")
            write_new_utils(slurm_utils_path, new)
            print(f"Updated {slurm_utils_path} for a cluster with no SMT.")


def main():
    args = get_args()
    print(f"Creating cluster profile in directory {args.output_dir}/slurm")

    # Download the profile, requires some user input
    get_cookiecutter_slurm_profile(output_dir=args.output_dir)

    # Check that parameters that need to be set a certain way are set that way
    check_profile_named_slurm(output_dir=args.output_dir)
    update_HoMi_reqs_in_config(output_dir=args.output_dir)
    check_use_conda_slurm(output_dir=args.output_dir)

    # Setup for cluster that doesn't have hyperthreading enabled
    if args.cluster_type=="slurm-nosmt":
        convert_slurm_profile_tasks(output_dir=args.output_dir)


if __name__=="__main__":
    main()