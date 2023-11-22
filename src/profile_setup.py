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
    print("Provide the following inputs to setup your cluster profile.\n")
    print("All default options can be used, EXCEPT:")
    print("****You MUST answer `True` to `use_conda`.****")
    print("****Please leave `profile_name` as `slurm`.**** \n") 

    cookiecutter("gh:Snakemake-Profiles/slurm", output_dir=output_dir)


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

    get_cookiecutter_slurm_profile(output_dir=args.output_dir)
    check_use_conda_slurm(output_dir=args.output_dir)

    if args.cluster_type=="slurm-nosmt":
        convert_slurm_profile_tasks(output_dir=args.output_dir)

if __name__=="__main__":
    main()