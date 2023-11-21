import argparse
import subprocess


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-type", 
                        choices=["slurm-smt", "slurm-nosmt"],
                        help="The type of cluster for which the profile should be created. \
                              Currently supported options are slurm-smt (slurm with simultaneous multithreading) \
                              and slurm-nosmt (slurm-managed cluster without simultaneous multithreading).",
                        required=True)
    return parser.parse_args()


def main():
    get_args()

if __name__=="__main__":
    main()