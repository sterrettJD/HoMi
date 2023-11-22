import argparse
import subprocess
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
                        required=True)
    return parser.parse_args()


def get_cookiecutter_slurm_profile(output_dir):
    print("Provide the following inputs to setup your cluster profile.")
    print("All default options can be used, EXCEPT:")
    print("****You MUST answer `True` to `use_conda`.**** \n") 
  
    cookiecutter("gh:Snakemake-Profiles/slurm", output_dir=output_dir)


def main():
    args = get_args()
    get_cookiecutter_slurm_profile(output_dir=args.output_dir)

if __name__=="__main__":
    main()