import argparse
import yaml
import pandas as pd
import os
from shutil import rmtree


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    parser.add_argument("metadata", help="the metadata file used for HoMi")
    parser.add_argument("--humann", action="store_true",
                        help="Pass this to clean up humann unneeded files")
    return parser.parse_args()


def read_config(config_filepath):
    with open(config_filepath, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


def read_sample_list(filepath):
    return pd.read_csv(filepath)["Sample"]


def clean_humann_temps(config, samples):
    prefix = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"
    humann_path = f"{prefix}.nonhost.humann"
    temp_dirs = [os.path.join(humann_path, sample, f"{sample}_humann_temp") for sample in samples]
    
    for i, sample in enumerate(samples):
        diamond_aligned = os.path.join(temp_dirs[i], f"{sample}_diamond_aligned.tsv")
        diamond_unaligned = os.path.join(temp_dirs[i], f"{sample}_diamond_unaligned.fa")
        print(f"Removing old file {diamond_aligned}")
        os.remove(diamond_aligned)
        print(f"Removing old file {diamond_unaligned}")
        os.remove(diamond_unaligned)

        remaining = os.listdir()
        other_temps = [x.startswith("tmp") for x in remaining]

        for dir in other_temps:
            print(f"Removing temporary directory {dir}")
            rmtree(dir)

    
def main():
    args = get_args()
    config = read_config(args.config)
    samples = read_sample_list(args.metadata)

    if args.humann is not None:
        clean_humann_temps(config, samples)



if __name__=="__main__":
    main()